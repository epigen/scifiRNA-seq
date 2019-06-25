#!/bin/bash

for i in "$@"
do
case $i in
    -n=*|--run-name=*)
    RUN_NAME="${i#*=}"
    shift # past argument=value
    ;;    --barcode_annotation=*)
    BARCODE_ANNOTATION="${i#*=}"
    shift # past argument=value
    ;;    --output-dir=*)
    ROOT_OUTPUT_DIR="${i#*=}"
    shift # past argument=value
    ;;    --cpus=*)
    CPUS="${i#*=}"
    shift # past argument=value
    ;;    --mem=*)
    MEM="${i#*=}"
    shift # past argument=value
    ;;    --queue=*)
    QUEUE="${i#*=}"
    shift # past argument=value
    ;;    --time=*)
    TIME="${i#*=}"
    shift # past argument=value
    ;;    --python-summarizer-script=*)
    PYTHON_SUMMARIZER_SCRIPT="${i#*=}"
    shift # past argument=value
    ;;    --round2_whitelist=*)
    ROUND2_WHITELIST="${i#*=}"
    shift # past argument=value
    ;;
    *)
          # unknown option
    ;;
esac
done

echo "RUN_NAME         = ${RUN_NAME}"
echo "ROOT DIRECTORY   = ${ROOT_OUTPUT_DIR}"
echo "SLURM PARAMETERS = $CPUS, $MEM, $QUEUE, $TIME"

# Start
mkdir -p $ROOT_OUTPUT_DIR
cd $ROOT_OUTPUT_DIR
echo "successfully entered directory"


# Summarize across lanes
echo "BARCODE_ANNOTATION  = ${BARCODE_ANNOTATION}"
for SAMPLE_NAME in `tail -n +2 $BARCODE_ANNOTATION | cut -d , -f 14`; do
JOB_NAME=scifi_pipeline.${SAMPLE_NAME}.summarize
SAMPLE_DIR=${ROOT_OUTPUT_DIR}/${SAMPLE_NAME}
JOB=${SAMPLE_DIR}/${JOB_NAME}.sh
LOG=${SAMPLE_DIR}/${JOB_NAME}.log

echo '#!/bin/env bash' > $JOB
echo "date" >> $JOB
echo "" >> $JOB
echo "python3 -u $PYTHON_SUMMARIZER_SCRIPT \
--sample-name $SAMPLE_NAME \
--r1-annot $BARCODE_ANNOTATION \
--r1-attributes plate plate_well donor_id sex \
--r2-barcodes $ROUND2_WHITELIST \
--cell-barcodes r2 \
--only-summary \
--no-save-intermediate \
--min-umi-output 20 \
--no-output-header \
--save-gene-expression \
${SAMPLE_DIR}/${SAMPLE_NAME}.*.STAR.Aligned.out.bam.featureCounts.bam \
${SAMPLE_DIR}/${SAMPLE_NAME}" >> $JOB
echo "" >> $JOB
echo "python3 -u $PYTHON_SUMMARIZER_SCRIPT \
--sample-name $SAMPLE_NAME \
--r1-annot $BARCODE_ANNOTATION \
--r1-attributes plate plate_well donor_id sex \
--r2-barcodes $ROUND2_WHITELIST \
--cell-barcodes r2 \
--only-summary \
--no-save-intermediate \
--min-umi-output 20 \
--no-output-header \
--save-gene-expression \
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
