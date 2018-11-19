
# Edit these parameters
RUN_NAMES=(
sci-RNA-seq_SCI012_12-nuclei_various-cycles
sci-RNA-seq_SCI012_24-nuclei_various-cycles
sci-RNA-seq_SCI012_12-nuclei_various-enzymes
sci-RNA-seq_SCI012_24-nuclei_various-enzymes
sci-RNA-seq_SCI013_12-cells_various-rtprimers
sci-RNA-seq_SCI013_12-nuclei_various-rtprimers
)
FLOWCELLS=(
BSF_0535_H2MJ2BGX9
BSF_0535_H2MJ2BGX9
BSF_0535_H2MJ2BGX9
BSF_0535_H2MJ2BGX9
BSF_0535_H2MJ2BGX9
BSF_0535_H2MJ2BGX9
)
BSF_NAMES=(
SCI_012_1_SSVI_12_nuclei
SCI_012_4_SSVI_24_nuclei
SCI_012_6_ProtoMax_12_nuclei
SCI_012_7_ProtoMax_24_nuclei
SCI_013_1_SSVI_12_cells
SCI_013_2_SSVI_12_nuclei
)

BARCODE_ANNOTATIONS=(
/scratch/lab_bock/shared/projects/sci-rna/metadata/sciRNA-seq.SCI012.oligos_2018-11-14.csv
/scratch/lab_bock/shared/projects/sci-rna/metadata/sciRNA-seq.SCI012.oligos_2018-11-14.csv
/scratch/lab_bock/shared/projects/sci-rna/metadata/sciRNA-seq.SCI012.oligos_2018-11-14.csv
/scratch/lab_bock/shared/projects/sci-rna/metadata/sciRNA-seq.SCI012.oligos_2018-11-14.csv
/scratch/lab_bock/shared/projects/sci-rna/metadata/sciRNA-seq.SCI013.oligos_2018-11-14.csv
/scratch/lab_bock/shared/projects/sci-rna/metadata/sciRNA-seq.SCI013.oligos_2018-11-14.csv
)


ROOT_OUTPUT_DIR=/scratch/lab_bock/shared/projects/sci-rna
N_PARTS=4
STEP=5000000
MAX_MISMATCHES=3
STAR_DIR=/data/groups/lab_bock/shared/resources/genomes/hg38_mm10_transgenes_Tcrlibrary/
GTF_FILE=/data/groups/lab_bock/shared/resources/genomes/hg38_mm10_transgenes_Tcrlibrary/Homo_sapiens-Mus_musculus.Ensembl92.dna.primary_assembly.Tcr_lambda_spiked.gtf


# Don't edit from here
mkdir -p $ROOT_OUTPUT_DIR/{logs,seqs,fastqc,star,barcodes,expression}
SAMPLE_NUMBERS=(`seq 0 $((${#RUN_NAMES[@]} - 1))`)
PARTS=`seq 1 $N_PARTS`


# Read to FASTQ
for SAMPLE in ${SAMPLE_NUMBERS[@]}; do
FLOWCELL=${FLOWCELLS[$SAMPLE]}
BSF_NAME=${BSF_NAMES[$SAMPLE]}
RUN_NAME=${RUN_NAMES[$SAMPLE]}
echo $FLOWCELL $BSF_NAME $RUN_NAME
for I in ${PARTS[@]}; do
sbatch -J ${RUN_NAME}.to_fastq.part_${I} \
-o ${ROOT_OUTPUT_DIR}/logs/${RUN_NAME}.to_fastq.part_${I}.log \
-c 8 --mem 80000 -p shortq \
--wrap "bedtools bamtofastq \
-i /scratch/lab_bsf/samples/${FLOWCELL}/${FLOWCELL}_${I}_samples/${FLOWCELL}_${I}#${BSF_NAME}.bam \
-fq ${ROOT_OUTPUT_DIR}/seqs/${RUN_NAME}.part_${I}.fastq"
done
done

# Map Reads
for SAMPLE in ${SAMPLE_NUMBERS[@]}; do
RUN_NAME=${RUN_NAMES[$SAMPLE]}
echo $RUN_NAME
for I in ${PARTS[@]}; do
sbatch -J ${RUN_NAME}.STAR.part_${I} \
-o ${ROOT_OUTPUT_DIR}/logs/${RUN_NAME}.STAR.part_${I}.log \
-c 12 --mem 200000 -p shortq \
--wrap "STAR --runThreadN 12 --genomeDir \
$STAR_DIR \
--outSAMunmapped Within --readFilesIn \
${ROOT_OUTPUT_DIR}/seqs/${RUN_NAME}.part_${I}.fastq \
--outFileNamePrefix ${ROOT_OUTPUT_DIR}/star/${RUN_NAME}.STAR.part_${I}. \
--outSAMtype BAM Unsorted --clip3pAdapterSeq AAAAAA"
done
done


# sort mapped/unmapped STAR output by name
for SAMPLE in ${SAMPLE_NUMBERS[@]}; do
RUN_NAME=${RUN_NAMES[$SAMPLE]}
echo $RUN_NAME
for I in ${PARTS[@]}; do
sbatch -J ${RUN_NAME}.sort.part_${I} \
-o ${ROOT_OUTPUT_DIR}/logs/${RUN_NAME}.sort.part_${I}.log \
-c 12 --mem 40000 -p shortq --time 07:30:00 \
--wrap "sambamba sort -t 8 -n \
${ROOT_OUTPUT_DIR}/star/${RUN_NAME}.STAR.part_${I}.Aligned.out.bam"
done
done


# tag with gene
for SAMPLE in ${SAMPLE_NUMBERS[@]}; do
RUN_NAME=${RUN_NAMES[$SAMPLE]}
echo $RUN_NAME
for I in ${PARTS[@]}; do
sbatch -J ${RUN_NAME}.htseq-count.part_${I} \
-o ${ROOT_OUTPUT_DIR}/logs/${RUN_NAME}.htseq-count.part_${I}.log \
-c 12 --mem 40000 -p shortq \
--wrap "samtools view \
${ROOT_OUTPUT_DIR}/star/${RUN_NAME}.STAR.part_${I}.Aligned.out.sorted.bam | \
~/.local/bin/htseq-count -f sam -a 10 -t exon -i gene_id \
--secondary-alignments=ignore --supplementary-alignments=ignore --additional-attr=gene_name \
--samout=${ROOT_OUTPUT_DIR}/star/${RUN_NAME}.STAR.htseq-count.sam.part_${I}.sam \
- \
$GTF_FILE > \
${ROOT_OUTPUT_DIR}/star/${RUN_NAME}.STAR.quant.part_${I}.tsv"
done
done


for SAMPLE in ${SAMPLE_NUMBERS[@]}; do
RUN_NAME=${RUN_NAMES[$SAMPLE]}
echo $RUN_NAME
for I in ${PARTS[@]}; do
# extract only read name
sbatch -J ${RUN_NAME}.extract_read.part_${I} \
-o ${ROOT_OUTPUT_DIR}/logs/${RUN_NAME}.extract_read.part_${I}.log \
-c 1 --mem 20000 -p shortq \
--wrap "\
cut -f 1 \
${ROOT_OUTPUT_DIR}/star/${RUN_NAME}.STAR.htseq-count.sam.part_${I}.sam \
> ${ROOT_OUTPUT_DIR}/star/${RUN_NAME}.STAR.htseq-count.read.part_${I}.txt"
# extract only gene
sbatch -J ${RUN_NAME}.extract_gene.part_${I} \
-o ${ROOT_OUTPUT_DIR}/logs/${RUN_NAME}.extract_gene.part_${I}.log \
-c 1 --mem 20000 -p shortq \
--wrap "\
sed 's/^HWI.*XF:Z://g' \
${ROOT_OUTPUT_DIR}/star/${RUN_NAME}.STAR.htseq-count.sam.part_${I}.sam \
> ${ROOT_OUTPUT_DIR}/star/${RUN_NAME}.STAR.htseq-count.gene.part_${I}.txt"
done
done


# join read and gene name
for SAMPLE in ${SAMPLE_NUMBERS[@]}; do
RUN_NAME=${RUN_NAMES[$SAMPLE]}
echo $RUN_NAME
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
done


# Extract sequences for all files
## with varying mismatches levels.
### This is only required with running the script in 'slim' mode,
### otherwise fat mode contains all matches which can be selected later 
for SAMPLE in ${SAMPLE_NUMBERS[@]}; do
FLOWCELL=${FLOWCELLS[$SAMPLE]}
BSF_NAME=${BSF_NAMES[$SAMPLE]}
RUN_NAME=${RUN_NAMES[$SAMPLE]}
BARCODE_ANNOTATION=${BARCODE_ANNOTATIONS[$SAMPLE]}
echo $FLOWCELL $BSF_NAME $RUN_NAME
for I in ${PARTS[@]}; do
SIZE=`samtools view -c /scratch/lab_bsf/samples/${FLOWCELL}/${FLOWCELL}_${I}_samples/${FLOWCELL}_${I}#${BSF_NAME}.bam`
for START in `seq 0 $STEP ${SIZE}`; do
END=$((START + STEP))
for MISMATCHES in `seq 0 $MAX_MISMATCHES`; do
echo $I $START $END $MISMATCHES
sbatch -J ${RUN_NAME}.barcode_extract.part_${I}.${START}_${END}.mis${MISMATCHES} \
-o ${ROOT_OUTPUT_DIR}/logs/${RUN_NAME}.barcode_extract.part_${I}.${START}_${END}.mis${MISMATCHES}.log \
-c 1 --mem 20000 -p shortq --time 1:00:00 \
--wrap "python -u ${ROOT_OUTPUT_DIR}/src/scirnaseq.extract_barcodes.py \
--mode slim --max-mismatches $MISMATCHES \
--start ${START} --end ${END} \
-a ${BARCODE_ANNOTATION} \
-o ${ROOT_OUTPUT_DIR}/barcodes/${RUN_NAME}.part_${I}.barcodes.${START}_${END}.mis_${MISMATCHES}.csv.gz \
/scratch/lab_bsf/samples/${FLOWCELL}/${FLOWCELL}_${I}_samples/${FLOWCELL}_${I}#${BSF_NAME}.bam"
# sleep 1
done
done
done
done


for SAMPLE in ${SAMPLE_NUMBERS[@]}; do
RUN_NAME=${RUN_NAMES[$SAMPLE]}
BARCODE_ANNOTATION=${BARCODE_ANNOTATIONS[$SAMPLE]}
for MISMATCHES in `seq 1 $((MAX_MISMATCHES+1))`; do
echo $RUN_NAME $((MISMATCHES-1)) $MISMATCHES
sbatch -J ${RUN_NAME}.inspect_barcodes.mis${MISMATCHES} \
-o ${ROOT_OUTPUT_DIR}/logs/${RUN_NAME}.inspect_barcodes.mis${MISMATCHES}.log \
-c 1 --mem 80000 -p shortq --time 1:00:00 \
--wrap "python -u ${ROOT_OUTPUT_DIR}/src/scirnaseq.inspect_barcodes.py \
--min-mismatches $((MISMATCHES-1)) --max-mismatches $MISMATCHES \
--barcodes round1,round2 \
--annotation ${BARCODE_ANNOTATION} \
$RUN_NAME"
done
done



# STATS
FEATURES=(
__not_aligned
__no_feature
__alignment_not_unique
ENSG
ENSMUS
lambda
Cas9_blast_gene
CTRL0
Tcrlibrary)

MISMATCHES=0

for SAMPLE in ${RUN_NAMES[@]}; do
echo "Sample" $SAMPLE
FLOWCELL=${FLOWCELLS[$SAMPLE]}
BSF_NAME=${BSF_NAMES[$SAMPLE]}
RUN_NAME=${RUN_NAMES[$SAMPLE]}
echo $FLOWCELL $BSF_NAME $RUN_NAME
for I in ${PARTS[@]}; do
echo "Part" $I
L=`samtools view -c /scratch/lab_bsf/samples/${FLOWCELL}/${FLOWCELL}_${I}_samples/${FLOWCELL}_${I}#${BSF_NAME}.bam`
echo $SAMPLE $I $MISMATCHES "sequenced_reads" $L >> ${ROOT_OUTPUT_DIR}/stats/${SAMPLE}.txt  # sequenced_reads

for FEATURE in ${FEATURES[@]}; do
echo "Feature" $FEATURE
L=`cat star/${SAMPLE}.STAR.htseq-count.gene.part_${I}.txt | grep $FEATURE | wc -l`
echo $SAMPLE $I $MISMATCHES $FEATURE $L >> ${ROOT_OUTPUT_DIR}/stats/${SAMPLE}.txt
done

L=`sed 's/,/\t/g' barcodes/${SAMPLE}.0mis.barcode_umi_dups.count.csv | cut -f 4 | paste -sd+ | bc`
echo $SAMPLE $I $MISMATCHES "total_umis" $L >> ${ROOT_OUTPUT_DIR}/stats/${SAMPLE}.txt  # reads sequenced

L=`sed 's/,/\t/g' barcodes/${SAMPLE}.0mis.barcode_gene_umi_count.clean.csv | cut -f 4 | paste -sd+ | bc`
echo $SAMPLE $I $MISMATCHES "umis" $L >> ${ROOT_OUTPUT_DIR}/stats/${SAMPLE}.txt  # total UMIs

done
done
