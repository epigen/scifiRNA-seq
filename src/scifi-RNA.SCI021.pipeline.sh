# Edit these parameters
## Sample
SAMPLE_NAME=sci-RNA-seq_SCI021_125K
FLOWCELL=BSF_0596_H7H5MDRXX
BSF_NAME=SCI021_125K_S50680
N_LANES=2
# SAMPLE_NAME=sci-RNA-seq_SCI020_3uL_reseq_4K
# FLOWCELL=BSF_0591_H5WTTBGXB
# BSF_NAME=SCI020_2_3ul_S50522
# N_LANES=4
BARCODE_ANNOTATION=/scratch/lab_bock/shared/projects/sci-rna/metadata/sciRNA-seq.SCI017.oligos_2019-02-11.csv

CPUS=12
MEM=120000
QUEUE=shortq

## project
ROOT_OUTPUT_DIR=/scratch/lab_bock/shared/projects/sci-rna
STAR_DIR=/home/arendeiro/resources/genomes/hg38_mm10_transgenes_Tcrlibrary/indexed_STAR_2.7.0e
GTF_FILE=/data/groups/lab_bock/shared/resources/genomes/hg38_mm10_transgenes_Tcrlibrary/Homo_sapiens-Mus_musculus.Ensembl92.dna.primary_assembly.Tcr_lambda_spiked.gtf
WHITELIST_10X=/scratch/lab_bock/shared/projects/sci-rna/metadata/737K-cratac-v1.reverse_complement.txt
# STEP=5000000
# MAX_MISMATCHES=3


# Don't edit from here
SAMPLE_DIR=${ROOT_OUTPUT_DIR}/data/${SAMPLE_NAME}
mkdir -p $SAMPLE_DIR #/{logs,fastqc,mapped,barcodes,expression}
cd $SAMPLE_DIR
PARTS=`seq 1 $N_LANES`


# # unfortunatelly, even though STAR can output to stdout 
# # and featureCounts read from stdin, one cannot pipe them as 
# # featureCounts does not support detailed BAM output with stdin


for PART in ${PARTS[@]}; do
INPUT_BAM=/scratch/lab_bsf/samples/${FLOWCELL}/${FLOWCELL}_${PART}_samples/${FLOWCELL}_${PART}#${BSF_NAME}.bam

JOB_NAME=scifiRNA-seq.${SAMPLE_NAME}.${PART}.filter
JOB=${SAMPLE_DIR}/${JOB_NAME}.sh
LOG=${SAMPLE_DIR}/${JOB_NAME}.log

echo '#!/bin/env bash' > $JOB

echo "date" >> $JOB

# align with STAR >=2.7.0e
echo "" >> $JOB
echo "STAR \
--runThreadN $CPUS \
--genomeDir $STAR_DIR \
--clip3pAdapterSeq AAAAAA \
--outSAMprimaryFlag AllBestScore \
--outSAMattributes All \
--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 \
--outSAMunmapped Within \
--outFileNamePrefix ${SAMPLE_DIR}/${SAMPLE_NAME}.${PART}.STAR. \
--outSAMtype BAM Unsorted \
--readFilesType SAM SE \
--readFilesCommand samtools view -h \
--quantMode TranscriptomeSAM GeneCounts \
--readFilesIn $INPUT_BAM" >> $JOB

# count all reads overlapping a gene
echo "" >> $JOB
echo "featureCounts \
-T $CPUS \
-F GTF \
-t gene \
-g gene_id \
--extraAttributes gene_name \
-Q 30 \
-s 0 \
-R BAM \
-a $GTF_FILE \
-o ${SAMPLE_DIR}/${SAMPLE_NAME}.${PART}.STAR.featureCounts.quant_gene.tsv \
${SAMPLE_DIR}/${SAMPLE_NAME}.${PART}.STAR.Aligned.out.bam" >> $JOB
# -M \
# -Q 0 \
# --primary \

# Filtering out secondary alignments and reads not assigned to a gene
# and convert to CSV
echo "" >> $JOB
echo "sambamba view \
-t $CPUS \
-F \"([XS] == 'Assigned') and (not unmapped) and (mapping_quality >= 30) and (not secondary_alignment) and (not supplementary)\" \
${SAMPLE_DIR}/${SAMPLE_NAME}.${PART}.STAR.Aligned.out.bam.featureCounts.bam \
| awk 'BEGIN {FS=\"\\t\"; OFS=\",\"} {print \$1, \$20, \$21, \$26, \$29, \$4}' \
| sed -r 's/.*:(.*:.*:.*):(.*)#.*,r1/\1:\2,r1/g' \
| sed 's/\w\w:Z://g' \
| gzip \
> ${SAMPLE_DIR}/${SAMPLE_NAME}.${PART}.STAR.filtered2.csv.gz" >> $JOB


echo "" >> $JOB
echo "date" >> $JOB
echo "" >> $JOB

sbatch -J $JOB_NAME \
-o $LOG \
-c $CPUS --mem $MEM -p $QUEUE \
$JOB
done


# Stats and summaries
# this uses Miller (https://github.com/johnkerl/miller)

INPUT_BAM=/scratch/lab_bsf/samples/${FLOWCELL}/${FLOWCELL}_${PART}_samples/${FLOWCELL}_${PART}#${BSF_NAME}.bam
FEATURECOUNTSFILE=${SAMPLE_DIR}/${SAMPLE_NAME}.stats.input_bam_files.txt
rm -f $FEATURECOUNTSFILE
for BAM_FILE in ${SAMPLE_DIR}/${SAMPLE_NAME}.*.STAR.Aligned.out.bam.featureCounts.bam; do
echo $BAM_FILE >> $FEATURECOUNTSFILE
done
readarray -t FEATURECOUNTSFILES < $FEATURECOUNTSFILE
STATS_FILE_PREFIX=${SAMPLE_DIR}/${SAMPLE_NAME}.stats_per_cell

CPUS=6
MEM=120000
QUEUE=longq
TIME=8-08:00:00
# # Count all reads per cell (from original BAM)
# # drawback: this includes non all barcodes! (must be filtered by correct ones)
JOB_NAME=${SAMPLE_NAME}.stats.raw_reads
LOG=${SAMPLE_DIR}/${JOB_NAME}.log
sbatch -J $JOB_NAME -c $CPUS --mem $MEM -o $LOG -p $QUEUE --time $TIME \
--wrap "samtools merge -@ 2 -f -u -O SAM /dev/stdout ${FEATURECOUNTSFILES[@]} \
| samtools view \
| mlr --ifs tab --ocsv \
    cut -f 1,20,21 \
    then rename 1,read,20,r1,21,r2 \
    then uniq -g r1,r2 -c -o sequenced_reads \
    then sort -f r1,r2 \
| sed 's/\w\w:Z://g' \
| gzip \
> ${STATS_FILE_PREFIX}.raw_reads.csv.gz"


# # Count unmapped reads (unique read names, not alignments) per cell
JOB_NAME=${SAMPLE_NAME}.stats.unmapped_reads
LOG=${SAMPLE_DIR}/${JOB_NAME}.log
sbatch -J $JOB_NAME -c $CPUS --mem $MEM -o $LOG -p $QUEUE --time $TIME \
--wrap "samtools merge -@ 2 -f -u -O SAM /dev/stdout ${FEATURECOUNTSFILES[@]} \
| sambamba view -S -F \"(unmapped)\" /dev/stdin \
| mlr --ifs tab --ocsv \
    cut -f 1,17,18 \
    then rename 1,read,17,r1,18,r2 \
    then uniq -g r1,r2 -c -o unmapped_reads \
    then sort -f r1,r2 \
| sed 's/\w\w:Z://g' \
| gzip \
> ${STATS_FILE_PREFIX}.unmapped_reads.csv.gz"

# # Count mapped, but unnassigned reads (unique read names, not alignments) per cell
JOB_NAME=${SAMPLE_NAME}.stats.mapped_unnassigned_reads
LOG=${SAMPLE_DIR}/${JOB_NAME}.log
sbatch -J $JOB_NAME -c $CPUS --mem $MEM -o $LOG -p $QUEUE --time $TIME \
--wrap "samtools merge -@ 2 -f -u -O SAM /dev/stdout ${FEATURECOUNTSFILES[@]} \
| sambamba view -S -F \"((not unmapped) and ([XS] != 'Assigned'))\" /dev/stdin \
| mlr --ifs tab --ocsv \
    cut -f 1,20,21 \
    then rename 1,read,20,r1,21,r2 \
    then cut -f read,r1,r2 \
    then uniq -g r1,r2,read \
    then uniq -g r1,r2 -c -o mapped_reads \
    then sort -f r1,r2 \
| sed 's/\w\w:Z://g' \
| gzip \
> ${STATS_FILE_PREFIX}.mapped_unnassigned_reads.csv.gz"

# # Count filtered reads (unique read names, not alignments) per cell
JOB_NAME=${SAMPLE_NAME}.stats.filtered_reads
LOG=${SAMPLE_DIR}/${JOB_NAME}.log
sbatch -J $JOB_NAME -c $CPUS --mem $MEM -o $LOG -p $QUEUE --time $TIME \
--wrap "samtools merge -@ 2 -f -u -O SAM /dev/stdout ${FEATURECOUNTSFILES[@]} \
| sambamba view -S -F \"([XS] == 'Assigned')\" /dev/stdin \
| mlr --ifs tab --ocsv \
    cut -f 1,20,21 \
    then rename 1,read,20,r1,21,r2 \
    then uniq -g r1,r2,read \
    then uniq -g r1,r2 -c -o filtered_reads \
    then sort -f r1,r2 \
| sed 's/\w\w:Z://g' \
| gzip \
> ${STATS_FILE_PREFIX}.filtered_reads.csv.gz"

# # Get duplication rate per cell (based on UMI, gene and position, then group per cell and count number of )
JOB_NAME=${SAMPLE_NAME}.stats.duplication
LOG=${SAMPLE_DIR}/${JOB_NAME}.log
sbatch -J $JOB_NAME -c $CPUS --mem $MEM -o $LOG -p $QUEUE --time $TIME \
--wrap "samtools merge -@ 2 -f -u -O SAM /dev/stdout ${FEATURECOUNTSFILES[@]} \
| sambamba view -S -F \"([XS] == 'Assigned')\" /dev/stdin \
| mlr --ifs tab --ocsv \
    cut -f 20,21,25,28,4 \
    then rename 20,r1,21,r2,25,umi,28,gene,4,pos \
    then count-distinct -f r1,r2,umi,gene,pos \
    then sort -f r1,r2 \
| sed 's/\w\w:Z://g' \
| gzip \
> ${STATS_FILE_PREFIX}.duplication.csv.gz"

# # Count UMIS per cell
JOB_NAME=${SAMPLE_NAME}.stats.umis
LOG=${SAMPLE_DIR}/${JOB_NAME}.log
sbatch -J $JOB_NAME -c $CPUS --mem $MEM -o $LOG -p $QUEUE --time $TIME \
--wrap "samtools merge -@ 2 -f -u -O SAM /dev/stdout ${FEATURECOUNTSFILES[@]} \
| sambamba view -S -F \"([XS] == 'Assigned')\" /dev/stdin \
| mlr --ifs tab --ocsv \
    cut -f 20,21,25,28,4 \
    then rename 20,r1,21,r2,25,umi,28,gene,4,pos \
    then uniq -g r1,r2 -c -o umis \
    then sort -f r1,r2 \
| sed 's/\w\w:Z://g' \
| gzip \
> ${STATS_FILE_PREFIX}.umis.csv.gz"

# # Count Genes per cell
JOB_NAME=${SAMPLE_NAME}.stats.genes
LOG=${SAMPLE_DIR}/${JOB_NAME}.log
sbatch -J $JOB_NAME -c $CPUS --mem $MEM -o $LOG -p $QUEUE --time $TIME \
--wrap "samtools merge -@ 2 -f -u -O SAM /dev/stdout ${FEATURECOUNTSFILES[@]} \
| sambamba view -S -F \"([XS] == 'Assigned')\" /dev/stdin \
| mlr --ifs tab --ocsv \
    cut -f 20,21,28 \
    then rename 20,r1,21,r2,28,gene \
    then uniq -g r1,r2 -c -o genes \
    then sort -f r1,r2 \
| sed 's/\w\w:Z://g' \
| gzip \
> ${STATS_FILE_PREFIX}.genes.csv.gz"

# # Get expression matrix
JOB_NAME=${SAMPLE_NAME}.stats.expression_matrix
LOG=${SAMPLE_DIR}/${JOB_NAME}.log
sbatch -J $JOB_NAME -c $CPUS --mem $MEM -o $LOG -p $QUEUE --time $TIME \
--wrap "samtools merge -@ 2 -f -u -O SAM /dev/stdout ${FEATURECOUNTSFILES[@]} \
| sambamba view -S -F \"([XS] == 'Assigned')\" /dev/stdin \
| mlr --ifs tab --ocsv \
    cut -f 20,21,25,28,4 \
    then rename 20,r1,21,r2,25,umi,28,gene,4,pos \
    then uniq -g r1,r2,gene -c -o umis \
    then sort -f r1,r2 \
| sed 's/\w\w:Z://g' \
| gzip \
> ${STATS_FILE_PREFIX}.expression_matrix.csv.gz"
