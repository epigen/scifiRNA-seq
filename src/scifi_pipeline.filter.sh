#!/bin/bash

for i in "$@"
do
case $i in
    -n=*|--run-name=*)
    RUN_NAME="${i#*=}"
    shift # past argument=value
    ;;
    -a=*|--annotation=*)
    BARCODE_ANNOTATION="${i#*=}"
    shift # past argument=value
    ;;
    -v=*|--variables=*)
    VARIABLES="${i#*=}"
    shift # past argument=value
    ;;
    --species-mixture=*)
    SPECIES_MIXTURE="${i#*=}"
    shift # past argument=value
    ;;
    --output-dir=*)
    ROOT_OUTPUT_DIR="${i#*=}"
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
echo "BARCODE_ANNOTATION = ${BARCODE_ANNOTATION}"
echo "VARIABLES        = ${VARIABLES}"
echo "SPECIES_MIXTURE  = ${SPECIES_MIXTURE}"
echo "ROOT DIRECTORY   = ${ROOT_OUTPUT_DIR}"
echo "SLURM PARAMETERS = $CPUS, $MEM, $QUEUE, $TIME"

# Start

# # Make path absolute
if [[ ! "$BARCODE_ANNOTATION" = /* ]]; then
    BARCODE_ANNOTATION=`pwd`/${BARCODE_ANNOTATION}
fi

SUMMARIZER=`pwd`/src/scifi_pipeline.summarizer.py


ADDITIONAL_ARGS=""
if [[ $SPECIES_MIXTURE = "1" ]]; then
    ADDITIONAL_ARGS="$ADDITIONAL_ARGS --species-mixture "
fi

mkdir -p $ROOT_OUTPUT_DIR
cd $ROOT_OUTPUT_DIR


# Summarize across lanes
for SAMPLE_NAME in `tail -n +2 $BARCODE_ANNOTATION | cut -d , -f 1`; do
JOB_NAME=scifi_pipeline.${SAMPLE_NAME}.summarize
SAMPLE_DIR=${ROOT_OUTPUT_DIR}/${SAMPLE_NAME}
JOB=${SAMPLE_DIR}/${JOB_NAME}.sh
LOG=${SAMPLE_DIR}/${JOB_NAME}.log

echo '#!/bin/env bash' > $JOB
echo "date" >> $JOB
echo "" >> $JOB
echo "python3 -u $SUMMARIZER \
--sample-name $SAMPLE_NAME \
--r1-annot $BARCODE_ANNOTATION \
--r1-attributes $VARIABLES \
--cell-barcodes r2 \
--only-summary \
--no-save-intermediate \
--min-umi-output 20 \
--no-output-header \
--save-gene-expression \
$ADDITIONAL_ARGS \
${SAMPLE_DIR}/${SAMPLE_NAME}.*.STAR.Aligned.out.bam.featureCounts.bam \
${SAMPLE_DIR}/${SAMPLE_NAME}" >> $JOB
echo "" >> $JOB
echo "date" >> $JOB
echo "" >> $JOB

sbatch -J $JOB_NAME \
-o $LOG --time $TIME \
-c $CPUS --mem $MEM -p $QUEUE \
$JOB


JOB_NAME=scifi_pipeline.${SAMPLE_NAME}.summarize-exon
SAMPLE_DIR=${ROOT_OUTPUT_DIR}/${SAMPLE_NAME}
JOB=${SAMPLE_DIR}/${JOB_NAME}.sh
LOG=${SAMPLE_DIR}/${JOB_NAME}.log

echo '#!/bin/env bash' > $JOB
echo "date" >> $JOB
echo "" >> $JOB
echo "python3 -u $SUMMARIZER \
--sample-name $SAMPLE_NAME \
--r1-annot $BARCODE_ANNOTATION \
--r1-attributes $VARIABLES \
--cell-barcodes r2 \
--only-summary \
--no-save-intermediate \
--min-umi-output 20 \
--no-output-header \
--save-gene-expression \
$ADDITIONAL_ARGS \
${SAMPLE_DIR}/${SAMPLE_NAME}.*.STAR.Aligned.out.exon.bam.featureCounts.bam \
${SAMPLE_DIR}/${SAMPLE_NAME}.exon" >> $JOB
echo "" >> $JOB
echo "date" >> $JOB
echo "" >> $JOB

sbatch -J $JOB_NAME \
-o $LOG --time $TIME \
-c $CPUS --mem $MEM -p $QUEUE \
$JOB
done
