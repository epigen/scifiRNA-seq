

# Edit these parameters

ROOT_OUTPUT_DIR=/scratch/lab_bock/shared/projects/sci-rna/
RUN_NAME=BSF_0477_HJ7J7BGX2
N_PARTS=4
STEP=5000000
MAX_MISMATCHES=1

STAR_DIR=/data/groups/lab_bock/shared/resources/genomes/hg38_mm10_transgenes_Tcrlibrary/
GTF_FILE=/data/groups/lab_bock/shared/resources/genomes/hg38_mm10_transgenes_Tcrlibrary/Homo_sapiens-Mus_musculus.Ensembl92.dna.primary_assembly.Tcr_lambda_spiked.gtf
BARCODE_ANNOTATION=${ROOT_OUTPUT_DIR}/metadata/sciRNA-seq.oligos_2018-05-22v2.csv


# Don't edit from here

mkdir -p $ROOT_OUTPUT_DIR/{logs,seqs,fastqc,star,barcodes,expression}
PARTS=`seq 1 $N_PARTS`


# Barcodes tag to FASTQ
for I in ${PARTS[@]}; do
sbatch -J ${RUN_NAME}.to_fastq.part_${I} \
-o ${ROOT_OUTPUT_DIR}/logs/${RUN_NAME}.to_fastq.part_${I}.log \
-c 8 --mem 80000 -p shortq \
--wrap "samtools view \
/data/groups/lab_bsf/sequences/${RUN_NAME}/${RUN_NAME}_${I}.bam | \
cut -f 1,12,14 | \
grep -v 'RG:Z:' | \
grep -v 'XN' | \
sed 's/^/@/g' | \
sed 's/\t/\n/g' | \
sed 's/BC:Z://g' | \
sed 's/QT:Z:/+\n/g' | \
gzip > ${ROOT_OUTPUT_DIR}/seqs/${RUN_NAME}.part_${I}.fastq.gz"
done

# Fastqc only on barcode read
for I in ${PARTS[@]}; do
sbatch -J ${RUN_NAME}.fastqc.part_${I}.barcode_read_only \
-o ${ROOT_OUTPUT_DIR}/logs/${RUN_NAME}.fastqc.part_${I}.barcode_read_only.log \
-c 8 --mem 80000 -p shortq \
--wrap "fastqc -t 8 --noextract -o ${ROOT_OUTPUT_DIR}/fastqc \
${ROOT_OUTPUT_DIR}/seqs/${RUN_NAME}.part_${I}.fastq.gz"
done

# Read2 to FASTQ
for I in ${PARTS[@]}; do
sbatch -J ${RUN_NAME}.read2_to_bam.part_${I} \
-o ${ROOT_OUTPUT_DIR}/logs/${RUN_NAME}.read2_to_bam.part_${I}.log \
-c 8 --mem 80000 -p shortq \
--wrap "samtools view -@ 8 -f 141 -b -o \
${ROOT_OUTPUT_DIR}/seqs/${RUN_NAME}_${I}.read2.bam \
/data/groups/lab_bsf/sequences/${RUN_NAME}/${RUN_NAME}.part_${I}.bam"
done

# Fastqc only on read2
for I in ${PARTS[@]}; do
sbatch -J ${RUN_NAME}.fastqc.read2_only.part_${I} \
-o ${ROOT_OUTPUT_DIR}/logs/${RUN_NAME}.fastqc.read2_only.part_${I}.log \
-c 8 --mem 80000 -p shortq \
--wrap "fastqc -t 8 --noextract -o ${ROOT_OUTPUT_DIR}/fastqc \
${ROOT_OUTPUT_DIR}/seqs/${RUN_NAME}.read2.part_${I}.bam"
done


# Read1 to FASTQ
for I in ${PARTS[@]}; do
sbatch -J ${RUN_NAME}.read1_to_bam.part_${I} \
-o ${ROOT_OUTPUT_DIR}/logs/${RUN_NAME}.read1_to_bam.part_${I}.log \
-c 8 --mem 80000 -p shortq \
--wrap "samtools view -@ 8 -f 77 -b -o \
${ROOT_OUTPUT_DIR}/seqs/${RUN_NAME}.part_${I}.read1.bam \
/data/groups/lab_bsf/sequences/${RUN_NAME}/${RUN_NAME}.part_${I}.bam"
done


# Fastqc only on read1
for I in ${PARTS[@]}; do
sbatch -J ${RUN_NAME}.fastqc.read2_only.part_${I} \
-o ${ROOT_OUTPUT_DIR}/logs/${RUN_NAME}.fastqc.read2_only.part_${I}.log \
-c 8 --mem 80000 -p shortq \
--wrap "fastqc -t 8 --noextract \
-o ${ROOT_OUTPUT_DIR}/fastqc \
${ROOT_OUTPUT_DIR}/seqs/${RUN_NAME}.part_${I}.read1.bam"
done

# Map read2
for I in ${PARTS[@]}; do
sbatch -J ${RUN_NAME}.read2_to_fastq.part_${I} \
-o ${ROOT_OUTPUT_DIR}/logs/${RUN_NAME}.read2_to_fastq.part_${I}.log \
-c 6 --mem 80000 -p shortq \
--wrap "samtools view -f 141 \
${ROOT_OUTPUT_DIR}/seqs/${RUN_NAME}.read2.part_${I}.bam | \
cut -f 1,10,11 | \
sed 's/^/@/g' | \
sed 's/\t/\n/1' | \
sed 's/\t/\n+\n/1' | \
gzip > ${ROOT_OUTPUT_DIR}/seqs/${RUN_NAME}.read2.part_${I}.fastq.gz"
done

for I in ${PARTS[@]}; do
sbatch -J ${RUN_NAME}.read2_to_fastq.part_${I} \
-o ${ROOT_OUTPUT_DIR}/logs/${RUN_NAME}.read2_to_fastq.part_${I}.log \
-c 1 --mem 80000 -p shortq \
--wrap "bedtools bamtofastq -i \
${ROOT_OUTPUT_DIR}/seqs/${RUN_NAME}.read2.part_${I}.bam \
-fq ${ROOT_OUTPUT_DIR}/seqs/${RUN_NAME}.read2.part_${I}.fastq"
done

for I in ${PARTS[@]}; do
sbatch -J ${RUN_NAME}.STAR.part_${I} \
-o ${ROOT_OUTPUT_DIR}/logs/${RUN_NAME}.STAR.part_${I}.log \
-c 12 --mem 200000 -p longq \
--wrap "STAR --runThreadN 12 --genomeDir \
$STAR_DIR \
--outSAMunmapped Within --readFilesIn \
${ROOT_OUTPUT_DIR}/seqs/${RUN_NAME}.read2.part_${I}.fastq \
--outFileNamePrefix ${ROOT_OUTPUT_DIR}/star/${RUN_NAME}.STAR.part_${I}. \
--outSAMtype BAM Unsorted --clip3pAdapterSeq AAAAAA"
done

# sort mapped/unmapped STAR output by name
for I in ${PARTS[@]}; do
sbatch -J ${RUN_NAME}.sort.part_${I} \
-o ${ROOT_OUTPUT_DIR}/logs/${RUN_NAME}.sort.part_${I}.log \
-c 12 --mem 40000 -p longq \
--wrap "sambamba sort -t 8 -n \
${ROOT_OUTPUT_DIR}/star/${RUN_NAME}.STAR.aligned.part_${I}.bam"
done


# tag with gene
for I in ${PARTS[@]}; do
sbatch -J ${RUN_NAME}.htseq-count.part_${I} \
-o ${ROOT_OUTPUT_DIR}/logs/${RUN_NAME}.htseq-count.part_${I}.log \
-c 12 --mem 40000 -p longq \
--wrap "samtools view \
${ROOT_OUTPUT_DIR}/star/${RUN_NAME}.STAR.aligned.part_${I}.sorted.bam | \
htseq-count -f sam -a 10 -t exon -i gene_id \
--secondary-alignments=ignore --supplementary-alignments=ignore --additional-attr=gene_name \
--samout=${ROOT_OUTPUT_DIR}/star/${RUN_NAME}.STAR.htseq-count.sam -.part_${I}.sam \
$GTF_FILE > \
${ROOT_OUTPUT_DIR}/star/${RUN_NAME}.STAR.quant.part_${I}.tsv"
done

# extract only read name
for I in ${PARTS[@]}; do
sbatch -J ${RUN_NAME}.extract_read.part_${I} \
-o ${ROOT_OUTPUT_DIR}/logs/${RUN_NAME}.extract_read.part_${I}.log \
-c 1 --mem 20000 -p shortq \
--wrap "\
cut -f 1 \
${ROOT_OUTPUT_DIR}/star/${RUN_NAME}.STAR.htseq-count.sam.part_${I}.sam \
> ${ROOT_OUTPUT_DIR}/star/${RUN_NAME}.STAR.htseq-count.read.part_${I}.txt"
done
# extract only gene
for I in ${PARTS[@]}; do
sbatch -J ${RUN_NAME}.extract_gene.part_${I} \
-o ${ROOT_OUTPUT_DIR}/logs/${RUN_NAME}.extract_gene.part_${I}.log \
-c 1 --mem 20000 -p shortq \
--wrap "\
sed 's/^HWI.*XF:Z://g' \
${ROOT_OUTPUT_DIR}/star/${RUN_NAME}.STAR.htseq-count.sam.part_${I}.sam \
> ${ROOT_OUTPUT_DIR}/star/${RUN_NAME}.STAR.htseq-count.gene.part_${I}.txt"
done
# join read and gene name
for I in ${PARTS[@]}; do
sbatch -J ${RUN_NAME}.join_read_gene.part_${I} \
-o ${ROOT_OUTPUT_DIR}/logs/${RUN_NAME}.join_read_gene.part_${I}.log \
-c 1 --mem 20000 -p shortq \
--wrap "\
paste \
${ROOT_OUTPUT_DIR}/star/${RUN_NAME}.STAR.htseq-count.read.part_${I}.txt \
${ROOT_OUTPUT_DIR}/star/${RUN_NAME}.STAR.htseq-count.gene.part_${I}.txt \
| sed 's/\t/,/g' \
> ${ROOT_OUTPUT_DIR}/star/${RUN_NAME}.STAR.htseq-count.read_gene.part_${I}.csv"
done


# Extract sequences for all files
## with varying mismatches levels.
### This is only required with running the script in 'slim' mode,
### otherwise fat mode contains all matches which can be selected later 
for MISMATCHES in `seq 0 $MAX_MISMATCHES`; do
for I in ${PARTS[@]}; do
SIZE=`samtools view -c ${ROOT_OUTPUT_DIR}/seqs/${RUN_NAME}_${I}.read1.bam`
for START in `seq 0 $STEP ${SIZE}`; do
END=$((START + STEP))
echo $I $START $END
sbatch -J ${RUN_NAME}.barcode_extract.part_${I}.${START}_${END}.mis1 \
-o ${ROOT_OUTPUT_DIR}/logs/${RUN_NAME}.barcode_extract.part_${I}.${START}_${END}.mis1.log \
-c 1 --mem 20000 -p shortq --time 1:00:00 \
--wrap "python -u ~/scirnaseq.barcodes.py \
--mode slim --max-mismatches $MISMATCHES \
--start ${START} --end ${END} \
-a ${BARCODE_ANNOTATION} \
-o ${ROOT_OUTPUT_DIR}/barcodes/${RUN_NAME}.part_${I}.barcodes.${START}_${END}.mis_${MISMATCHES}.csv.gz \
${ROOT_OUTPUT_DIR}/seqs/${RUN_NAME}.part_${I}.read1.bam"
# sleep 1
done
done
done


# Run inspect_barcodes.py
## inputs:
## - ${ROOT_OUTPUT_DIR}/star/${RUN_NAME}.STAR.htseq-count.read_gene.part_${I}.csv
## - ${ROOT_OUTPUT_DIR}/barcodes/${RUN_NAME}.part_*.barcodes.*_*.mis_${MISMATCHES}.csv.gz
