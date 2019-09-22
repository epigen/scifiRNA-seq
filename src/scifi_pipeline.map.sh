#!/bin/bash

for i in "$@"
do
case $i in
    -n=*|--run-name=*)
    RUN_NAME="${i#*=}"
    shift # past argument=value
    ;;
    -f=*|--flowcell=*)
    FLOWCELL="${i#*=}"
    shift # past argument=value
    ;;
    -l=*|--n-lanes=*)
    N_LANES="${i#*=}"
    shift # past argument=value
    ;;
    -l=*|--n-barcodes=*)
    N_BARCODES="${i#*=}"
    shift # past argument=value
    ;;
    -a=*|--annotation=*)
    BARCODE_ANNOTATION="${i#*=}"
    shift # past argument=value
    ;;
    --output-dir=*)
    ROOT_OUTPUT_DIR="${i#*=}"
    shift # past argument=value
    ;;
    --star-exe=*)
    STAR_EXE="${i#*=}"
    shift # past argument=value
    ;;
    --star-dir=*)
    STAR_DIR="${i#*=}"
    shift # past argument=value
    ;;
    --gtf=*)
    GTF_FILE="${i#*=}"
    shift # past argument=value
    ;;
    --cpus=*)
    CPUS="${i#*=}"
    shift # past argument=value
    ;;
    --mem=*)
    MEM="${i#*=}"
    shift # past argument=value
    ;;
    --queue=*)
    QUEUE="${i#*=}"
    shift # past argument=value
    ;;
    --time=*)
    TIME="${i#*=}"
    shift # past argument=value
    ;;
    *)
          # unknown option
    ;;
esac
done

echo "RUN_NAME         = ${RUN_NAME}"
echo "FLOWCELL         = ${FLOWCELL}"
echo "NUMBER OF LANES  = ${N_LANES}"
echo "BARCODE NUMBER   = ${N_BARCODES}"
echo "ANNOTATION       = ${BARCODE_ANNOTATION}"
echo "ROOT DIRECTORY   = ${ROOT_OUTPUT_DIR}"
echo "STAR EXECUTABLE  = ${STAR_EXE}"
echo "STAR DIRECTORY   = ${STAR_DIR}"
echo "GTF FILE         = ${GTF_FILE}"
echo "SLURM PARAMETERS = $CPUS, $MEM, $QUEUE, $TIME"

# Don't edit from here

if [[ ! "$BARCODE_ANNOTATION" = /* ]]; then
    BARCODE_ANNOTATION=`pwd`/${BARCODE_ANNOTATION}
fi

mkdir -p $ROOT_OUTPUT_DIR
cd $ROOT_OUTPUT_DIR
LANES=`seq 1 $N_LANES`


# # unfortunatelly, even though STAR can output to stdout 
# # and featureCounts read from stdin, one cannot pipe them as 
# # featureCounts does not support detailed BAM output with stdin
for SAMPLE_NAME in `tail -n +2 $BARCODE_ANNOTATION | cut -d , -f 1`; do

# for LANE in ${LANES[@]}; do
SAMPLE_DIR=${ROOT_OUTPUT_DIR}/${SAMPLE_NAME}
mkdir -p $SAMPLE_DIR #/{logs,fastqc,mapped,barcodes,expression}

function join_by { local IFS="$1"; shift; echo "$*"; }

if [[ $string == *"gRNA"* ]]; then
    INPUT_BAM=${SAMPLE_DIR}/${SAMPLE_NAME}.${LANE}.trimmed.bam
else
    BAMS=()
    DIR=/scratch/users/dbarreca/private/custom_demux/scRNA/${FLOWCELL}
    # DIR=/scratch/users/arendeiro/custom_demux/scRNA
    # DIR=/scratch/lab_bsf/samples/${FLOWCELL}
    if [[ $N_BARCODES -gt 1 ]]; then
        for LANE in ${LANES[@]}; do
            BAMS+=(`eval echo $DIR/${FLOWCELL}_${LANE}_samples/${FLOWCELL}_${LANE}#${SAMPLE_NAME}_{01..$N_BARCODES}.bam`)
        done
    else
        for LANE in ${LANES[@]}; do
            BAMS+=(`echo $DIR/${FLOWCELL}_${LANE}_samples/${FLOWCELL}_${LANE}#${SAMPLE_NAME}.bam`)
        done
    fi
    LANE="ALL"
fi

INPUT_BAM=`join_by , ${BAMS[@]}`


JOB_NAME=scifi_pipeline.${SAMPLE_NAME}.${LANE}.map
JOB=${SAMPLE_DIR}/${JOB_NAME}.sh
LOG=${SAMPLE_DIR}/${JOB_NAME}.log

echo '#!/bin/env bash' > $JOB

echo "date" >> $JOB
echo "" >> $JOB

echo "#RUN_NAME         = ${RUN_NAME}" >> $JOB
echo "#FLOWCELL         = ${FLOWCELL}" >> $JOB
echo "#NUMBER OF LANES  = ${N_LANES}" >> $JOB
echo "#BARCODE NUMBER   = ${N_BARCODES}" >> $JOB
echo "#ANNOTATION       = ${BARCODE_ANNOTATION}" >> $JOB
echo "#ROOT DIRECTORY   = ${ROOT_OUTPUT_DIR}" >> $JOB
echo "#STAR EXECUTABLE  = ${STAR_EXE}" >> $JOB
echo "#STAR DIRECTORY   = ${STAR_DIR}" >> $JOB
echo "#GTF FILE         = ${GTF_FILE}" >> $JOB
echo "#SLURM PARAMETERS = $CPUS, $MEM, $QUEUE, $TIME" >> $JOB


# align with STAR >=2.7.0e
echo "" >> $JOB
echo "$STAR_EXE \
--runThreadN $CPUS \
--genomeDir $STAR_DIR \
--clip3pAdapterSeq AAAAAA \
--outSAMprimaryFlag AllBestScore \
--outSAMattributes All \
--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 \
--outSAMunmapped Within \
--outFileNamePrefix ${SAMPLE_DIR}/${SAMPLE_NAME}.${LANE}.STAR. \
--outSAMtype BAM Unsorted \
--readFilesType SAM SE \
--readFilesCommand samtools view -h \
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
-o ${SAMPLE_DIR}/${SAMPLE_NAME}.${LANE}.STAR.featureCounts.quant_gene.tsv \
${SAMPLE_DIR}/${SAMPLE_NAME}.${LANE}.STAR.Aligned.out.bam" >> $JOB

# Same as above but just for exons
echo "" >> $JOB
echo "ln -s ${SAMPLE_DIR}/${SAMPLE_NAME}.${LANE}.STAR.Aligned.out.bam \
${SAMPLE_DIR}/${SAMPLE_NAME}.${LANE}.STAR.Aligned.out.exon.bam" >> $JOB
echo "featureCounts \
-T $CPUS \
-F GTF \
-t exon \
-g gene_id \
--extraAttributes gene_name \
-Q 30 \
-s 0 \
-R BAM \
-a $GTF_FILE \
-o ${SAMPLE_DIR}/${SAMPLE_NAME}.${LANE}.STAR.featureCounts.quant_gene.exon.tsv \
${SAMPLE_DIR}/${SAMPLE_NAME}.${LANE}.STAR.Aligned.out.exon.bam" >> $JOB

echo "" >> $JOB
echo "date" >> $JOB
echo "" >> $JOB

sbatch -J $JOB_NAME \
-o $LOG --time $TIME \
-c $CPUS --mem $MEM -p $QUEUE \
$JOB

# done
done
