# Edit these parameters
## Sample
SAMPLE_NAME=SCI022_PBMC
FLOWCELL=BSF_0604_H5GLYBBXY
BSF_NAME=SCI022_PBMC
N_LANES=2
N_PARTS=4

SAMPLE_NAME=SCI022_TCell
FLOWCELL=BSF_0604_H5GLYBBXY
BSF_NAME=SCI022_TCell
N_LANES=2
N_PARTS=4

BARCODE_ANNOTATION=/scratch/lab_bock/shared/projects/sci-rna/metadata/sciRNA-seq.SCI022.oligos.csv

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
LANES=`seq 1 $N_LANES`
PARTS=`seq 1 $N_PARTS`


# # unfortunatelly, even though STAR can output to stdout 
# # and featureCounts read from stdin, one cannot pipe them as 
# # featureCounts does not support detailed BAM output with stdin

for LANE in ${LANES[@]}; do
for PART in ${PARTS[@]}; do
INPUT_BAM=/scratch/users/dbarreca/private/custom_demux/scRNA/${FLOWCELL}/${FLOWCELL}_${LANE}_samples/${FLOWCELL}_${LANE}#${BSF_NAME}_${PART}.bam

JOB_NAME=scifiRNA-seq.${SAMPLE_NAME}.${LANE}.${PART}
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
--outFileNamePrefix ${SAMPLE_DIR}/${SAMPLE_NAME}.${LANE}.${PART}.STAR. \
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
-o ${SAMPLE_DIR}/${SAMPLE_NAME}.${LANE}.${PART}.STAR.featureCounts.quant_gene.tsv \
${SAMPLE_DIR}/${SAMPLE_NAME}.${LANE}.${PART}.STAR.Aligned.out.bam" >> $JOB

# Filtering out secondary alignments and reads not assigned to a gene
# and convert to CSV
echo "" >> $JOB
echo "sambamba view \
-t $CPUS \
-F \"([XS] == 'Assigned') and (not unmapped) and (mapping_quality >= 30) and (not secondary_alignment) and (not supplementary)\" \
${SAMPLE_DIR}/${SAMPLE_NAME}.${LANE}.${PART}.STAR.Aligned.out.bam.featureCounts.bam \
| awk 'BEGIN {FS=\"\\t\"; OFS=\",\"} {print \$1, \$20, \$21, \$26, \$29, \$4}' \
| sed -r 's/.*:(.*:.*:.*):(.*)#.*,r1/\1:\2,r1/g' \
| sed 's/\w\w:Z://g' \
| gzip \
> ${SAMPLE_DIR}/${SAMPLE_NAME}.${LANE}.${PART}.STAR.filtered2.csv.gz" >> $JOB

echo "" >> $JOB
echo "date" >> $JOB
echo "" >> $JOB

sbatch -J $JOB_NAME \
-o $LOG \
-c $CPUS --mem $MEM -p $QUEUE \
$JOB
done
done
