# Preliminary analysis for sci-RNA-seq trial


# Parameters

SAMPLE_ANNOTATION=metadata/annotation.csv
CPUS=8
MEM=80000

# Vars
mkdir -p data
mkdir -p results

SAMPLE_NAMES=()
BAM_FILES=()

# Read annotation
counter=0
while IFS=, read -r col{1..17}
do
    (( counter++ ))
    [ $counter -eq 1 ] && continue  # skip first CSV line

    SAMPLE_NAMES+=($col1)
    BAM="/scratch/lab_bsf/samples/${col11}/${col11}_${col12}_samples/${col11}_${col12}#${col13}.bam"
    BAM_FILES+=($BAM)

    echo "Sample: '$col1'"
    echo "BAM file: '$BAM'"
    echo
done < $SAMPLE_ANNOTATION


# Pipeline
## Step 0: FASTQC & Flagstat
for ((I=0;I<${#SAMPLE_NAMES[@]};++i)); do
    mkdir -p data/${SAMPLE_NAMES[$I]}
    sambamba flagstat -t $CPUS ${BAM_FILES[$I]} > data/${SAMPLE_NAMES[$I]}/${SAMPLE_NAMES[$I]}.flagstat
    fastqc -t $CPUS -o data/${SAMPLE_NAMES[$I]} --no-extract ${BAM_FILES[$I]}
done


## Step 1: Trim reads, align, quantify transcripts
for ((I=0;I<${#SAMPLE_NAMES[@]};++i)); do
    mkdir -p data/${SAMPLE_NAMES[$I]}
    mkdir -p data/${SAMPLE_NAMES[$I]}/raw
    mkdir -p data/${SAMPLE_NAMES[$I]}/mapped

    # convert to FASTQ
    java -jar /cm/shared/apps/picard-tools/2.9.0/picard.jar SamToFastq \
    INPUT=${BAM_FILES[$I]} \
    FASTQ=data/${SAMPLE_NAMES[$I]}/raw/${SAMPLE_NAMES[$I]}.r1.fastq \
    SECOND_END_FASTQ=data/${SAMPLE_NAMES[$I]}/raw/${SAMPLE_NAMES[$I]}.r2.fastq \
    UNPAIRED_FASTQ=data/${SAMPLE_NAMES[$I]}/raw/${SAMPLE_NAMES[$I]}.unpaired.fastq

    # trim
    java -jar /cm/shared/apps/trimmomatic/0.36/trimmomatic-0.36.jar \
    SE -threads $CPUS -phred33 \
    data/${SAMPLE_NAMES[$I]}/raw/${SAMPLE_NAMES[$I]}.r2.fastq \
    data/${SAMPLE_NAMES[$I]}/raw/${SAMPLE_NAMES[$I]}.r2.trimmed.fastq \
    ILLUMINACLIP:/home/arendeiro/resources/adapters/illumina.fa:2:30:10 \
    LEADING:10 TRAILING:10 \
    SLIDINGWINDOW:4:15 MINLEN:32
    # CROP:50 \ # <- for comparisons with Quant-seq

    # align
    srun --mem $MEM -p develop /cm/shared/apps/star/2.5.2b/STAR \
    --outSAMtype BAM Unsorted \
    --runThreadN $CPUS \
    --genomeDir ~/projects/archive/crop-seq/spiked_genomes/hg38_spiked_Tcrlibrary/ \
    --readFilesIn data/${SAMPLE_NAMES[$I]}/raw/${SAMPLE_NAMES[$I]}.r2.trimmed.fastq \
    --outFileNamePrefix data/${SAMPLE_NAMES[$I]}/mapped/${SAMPLE_NAMES[$I]}.r2.trimmed.star \
    --outReadsUnmapped data/${SAMPLE_NAMES[$I]}/mapped/${SAMPLE_NAMES[$I]}.r2.trimmed.star.unmapped

    # sort, index BAM file
    sambamba sort -t $CPUS data/${SAMPLE_NAMES[$I]}/mapped/${SAMPLE_NAMES[$I]}.r2.trimmed.starAligned.out.bam

    # quantify transcripts
    htseq-count \
    -f bam -t gene -i gene_name \
    data/${SAMPLE_NAMES[$I]}/mapped/${SAMPLE_NAMES[$I]}.r2.trimmed.starAligned.out.sorted.bam \
    ~/projects/archive/crop-seq/spiked_genomes/hg38_spiked_Tcrlibrary/Homo_sapiens.GRCh38.77.spiked.gene_only.gtf \
    > data/${SAMPLE_NAMES[$I]}/${SAMPLE_NAMES[$I]}.r2.trimmed.starAligned.out.sorted.quant.tsv
done


# # Plot transcript coverage
# gtfToGenePred ~/projects/archive/crop-seq/spiked_genomes/hg38_spiked_Tcrlibrary/Homo_sapiens.GRCh38.77.spiked.gtf q
# awk '{print $2,$4,$5,$1,$8,$3}' q | sed 's/ /\t/g' | sortBed > ~/projects/archive/crop-seq/spiked_genomes/hg38_spiked_Tcrlibrary/Homo_sapiens.GRCh38.77.spiked.bed
for ((I=0;I<${#SAMPLE_NAMES[@]};++i)); do
    echo ${SAMPLE_NAMES[$I]}
    geneBody_coverage.py \
    -r ~/hg38_RefSeq_curated.ensembl_chroms.bed \
    -i data/${SAMPLE_NAMES[$I]}/mapped/${SAMPLE_NAMES[$I]}.r2.trimmed.starAligned.out.sorted.bam \
    -o genebody_coverage_${SAMPLE_NAMES[$I]}
done



# Plot replicate similarity
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv(os.path.join("metadata", "annotation.csv"))

quant = pd.DataFrame()
for sample in df['sample_name']:
    
    quant[sample] = pd.read_csv(os.path.join(
            "data", sample, sample + ".r2.trimmed.starAligned.out.sorted.quant.tsv"),
        sep="\t", squeeze=True, index_col=0, header=None)
quant = quant.loc[~quant.index.str.startswith("__")]
quant.to_csv(os.path.join("results", "quantification.raw_counts.csv"), index=True)

# normalize
quant = np.log2(1 + quant)
quant = (quant / quant.sum(axis=0)) * 1e5
quant.to_csv(os.path.join("results", "quantification.log2_rpm.csv"), index=True)


fig, axis = plt.subplots(1, 3, figsize=(3 * 4, 4), sharex=True, sharey=True)
axis[0].scatter(quant.iloc[:, 0], quant.iloc[:, 1], s=2, alpha=0.1, rasterized=True)
axis[1].scatter(quant.iloc[:, 0], quant.iloc[:, 2], s=2, alpha=0.1, rasterized=True)
axis[2].scatter(quant.iloc[:, 1], quant.iloc[:, 2], s=2, alpha=0.1, rasterized=True)
axis[0].plot([0, quant.max().max()], [0, quant.max().max()], "--", color="black", alpha=0.5, zorder=10)
axis[1].plot([0, quant.max().max()], [0, quant.max().max()], "--", color="black", alpha=0.5, zorder=10)
axis[2].plot([0, quant.max().max()], [0, quant.max().max()], "--", color="black", alpha=0.5, zorder=10)
axis[0].set_xlabel(quant.columns[0])
axis[0].set_ylabel(quant.columns[1])
axis[1].set_xlabel(quant.columns[0])
axis[1].set_ylabel(quant.columns[2])
axis[2].set_xlabel(quant.columns[1])
axis[2].set_ylabel(quant.columns[2])
sns.despine(fig)
fig.savefig(os.path.join("results", "sample_correlation.scatter.svg"), bbox_inches="tight", dpi=300)

fig, axis = plt.subplots(1, 3, figsize=(3 * 4, 4), sharex=False, sharey=False)
axis[0].scatter(quant.iloc[:, 0], quant.iloc[:, 1], s=2, alpha=0.1, rasterized=True)
axis[1].scatter(quant.iloc[:, 0], quant.iloc[:, 2], s=2, alpha=0.1, rasterized=True)
axis[2].scatter(quant.iloc[:, 1], quant.iloc[:, 2], s=2, alpha=0.1, rasterized=True)
axis[0].set_xlabel(quant.columns[0])
axis[0].set_ylabel(quant.columns[1])
axis[1].set_xlabel(quant.columns[0])
axis[1].set_ylabel(quant.columns[2])
axis[2].set_xlabel(quant.columns[1])
axis[2].set_ylabel(quant.columns[2])
sns.despine(fig)
fig.savefig(os.path.join("results", "sample_correlation.scatter.no_xy.svg"), bbox_inches="tight", dpi=300)
