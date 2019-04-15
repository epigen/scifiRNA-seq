#!/usr/bin/env python

import os


def write_lambda(
        output_fasta, output_gtf,
        seq_url=None,
        seq_name="lambda"):

    if seq_url is None:
        seq_url = "https://www.neb.com/-/media/nebus/page-images/tools-and-resources/interactive-tools/dna-sequences-and-maps/text-documents/lambdafsa.txt?la=en"
    fasta_header_template = ">{chrom}_chrom dna:chromosome chromosome:GRCh38:{chrom}_chrom:1:{length}:1 REF"

    gtf_template = """{chrom}_chrom\thavana\tgene\t1\t{length}\t.\t+\t.\tgene_id "{id}_gene"; gene_name "{id}_gene"; gene_source "ensembl_havana"; gene_biotype "lincRNA";
{chrom}_chrom\thavana\ttranscript\t1\t{length}\t.\t+\t.\tgene_id "{id}_gene"; transcript_id "{id}_transcript"; gene_name "{id}_gene"; gene_source "ensembl_havana"; gene_biotype "lincRNA"; transcript_name "{id}_transcript"; transcript_source "havana";
{chrom}_chrom\thavana\texon\t1\t{length}\t.\t+\t.\tgene_id "{id}_gene"; transcript_id "{id}_transcript"; exon_number "1"; gene_name "{id}_gene"; gene_source "ensembl_havana"; gene_biotype "lincRNA"; transcript_name "{id}_transcript"; transcript_source "havana"; exon_id "{id}_exon";
    """

    # get fasta
    os.system("curl {} > seq".format(seq_url))
    raw_sequence = open("seq", 'r').read().strip()
    fasta_sequence = ''.join(raw_sequence.split("\r\n")[1:]).replace("\r\n", "")
    fasta_header = fasta_header_template.format(chrom=seq_name, length=len(fasta_sequence))

    # get gtf
    gtf_entry = gtf_template.format(chrom=seq_name, id=seq_name, length=len(fasta_sequence))

    # write to file
    with open(output_fasta, "w") as fasta_handle:
        fasta_handle.writelines("\n".join([fasta_header, fasta_sequence]))
    with open(output_gtf, "w") as gtf_handle:
        gtf_handle.writelines(gtf_entry)


# initialize project
output_dir = "/home/arendeiro/resources/genomes"
genome_ref = "hg38_spiked"
library = "lambda"
spiked_dir = os.path.join(output_dir, genome_ref) + "_" + library

for _dir in [output_dir, spiked_dir]:
    if not os.path.exists(_dir):
        os.makedirs(_dir)

output_fasta = os.path.join(spiked_dir, "lambda.fa")
output_gtf = os.path.join(spiked_dir, "lambda.gtf")

# write gRNA library annotation
write_lambda(output_fasta, output_gtf)


# Make STAR index and supporting files
# Get fasta genome
os.system(
    "wget -O {} ftp://ftp.ensembl.org/pub/release-77/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    .format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz")))
os.system("gzip -d {}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz")))
# Get gtf annotation
os.system(
    "wget -O {} ftp://ftp.ensembl.org/pub/release-77/gtf/homo_sapiens/Homo_sapiens.GRCh38.77.gtf.gz"
    .format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.77.gtf.gz")))
os.system("gzip -d {}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.77.gtf.gz")))

# Add extra chromosomes (constructs) to genome
cmd = (
    "cat {} {} > {}"
    .format(
        os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.fa"),
        output_fasta,
        os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.spiked.fa")))
os.system(cmd)
cmd = (
    "cat {} {} > {}"
    .format(
        os.path.join(spiked_dir, "Homo_sapiens.GRCh38.77.gtf"),
        output_gtf,
        os.path.join(spiked_dir, "Homo_sapiens.GRCh38.77.spiked.gtf")))
os.system(cmd)

# Build STAR index (contruct + spiked with gRNAs)
cmd = "srun --mem 80000 -p develop /cm/shared/apps/star/2.4.2a/STAR"
cmd += " --runThreadN 8"
cmd += " --runMode genomeGenerate"
cmd += " --genomeDir {}".format(spiked_dir)
cmd += " --genomeFastaFiles {}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.spiked.fa"))
cmd += " --sjdbGTFfile {}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.77.spiked.gtf"))
cmd += " --sjdbOverhang 74"
os.system(cmd)

# Create sequence dictionaries (for piccard)
cmd = "srun --mem 80000 -p develop java -Xmx8g -jar /cm/shared/apps/picard-tools/1.140/picard.jar"
cmd += " CreateSequenceDictionary"
cmd += " REFERENCE={}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.spiked.fa"))
cmd += " OUTPUT={}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.spiked.dict"))
cmd += " GENOME_ASSEMBLY={}".format(genome_ref)
cmd += " SPECIES=human"
os.system(cmd)

# Create reflat files
cmd = "srun --mem 80000 -p develop java -Xmx80g -jar ~/Drop-seq_tools-1.12/jar/dropseq.jar ConvertToRefFlat"
cmd += " SEQUENCE_DICTIONARY={}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.spiked.dict"))
cmd += " ANNOTATIONS_FILE= {}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.77.spiked.gtf"))
cmd += " OUTPUT={}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.spiked.refFlat"))
os.system(cmd)

# Remove vanilla genome
os.remove(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.fa"))
os.remove(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.77.gtf"))


# Build STAR index for newer version (contruct + spiked with gRNAs)
cmd = "srun -J new_STAR_reference --mem 80000 -p develop /workspace/STAR-2.7.0e/bin/Linux_x86_64_static/STAR"
cmd += " --runThreadN 8"
cmd += " --runMode genomeGenerate"
cmd += " --genomeDir {}".format(spiked_dir)
cmd += " --genomeFastaFiles {}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.spiked.fa"))
cmd += " --sjdbGTFfile {}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.77.spiked.gtf"))
cmd += " --sjdbOverhang 74"
os.system(cmd)

# sbatch \
# -J new_STAR_reference \
# -c 12 --mem 80000 -p develop \
# --wrap "/cm/shared/apps/star/2.4.2a/STAR \
# --runThreadN 8 \
# --runMode genomeGenerate \
# --genomeDir /home/arendeiro/resources/genomes/hg38_spiked_lambda/indexed_star_index_2.4.2a \
# --genomeFastaFiles /home/arendeiro/resources/genomes/hg38_spiked_lambda/Homo_sapiens.GRCh38.dna.primary_assembly.spiked.fa \
# --sjdbGTFfile /home/arendeiro/resources/genomes/hg38_spiked_lambda/Homo_sapiens.GRCh38.77.spiked.gtf \
# --sjdbOverhang 74"


# sbatch \
# -J new_STAR_reference \
# -c 12 --mem 80000 -p develop \
# --wrap "~/workspace/STAR-2.7.0e/bin/Linux_x86_64_static/STAR \
# --runThreadN 12 \
# --runMode genomeGenerate \
# --genomeDir /home/arendeiro/resources/genomes/hg38_mm10_transgenes_Tcrlibrary/indexed_STAR_2.7.0e \
# --genomeFastaFiles /home/arendeiro/resources/genomes/hg38_mm10_transgenes_Tcrlibrary/Homo_sapiens-Mus_musculus.Ensembl92.dna.primary_assembly.Tcr_lambda_spiked.fa \
# --sjdbGTFfile /home/arendeiro/resources/genomes/hg38_mm10_transgenes_Tcrlibrary/Homo_sapiens-Mus_musculus.Ensembl92.dna.primary_assembly.Tcr_lambda_spiked.gtf \
# --sjdbOverhang 74"
