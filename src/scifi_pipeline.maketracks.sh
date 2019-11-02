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
    --array-size=*)
    ARRAY_SIZE="${i#*=}"
    shift # past argument=value
    ;;
    *)
          # unknown option
    ;;
esac
done

echo "RUN_NAME             = ${RUN_NAME}"
echo "BARCODE_ANNOTATION   = ${BARCODE_ANNOTATION}"
echo "ROOT DIRECTORY       = ${ROOT_OUTPUT_DIR}"
echo "SLURM PARAMETERS     = $CPUS, $MEM, $QUEUE, $TIME, $ARRAY_SIZE"

# Start

# # Make path absolute
if [[ ! "$BARCODE_ANNOTATION" = /* ]]; then
    BARCODE_ANNOTATION=`pwd`/${BARCODE_ANNOTATION}
fi

JOB_DESCR="maketracks"


mkdir -p $ROOT_OUTPUT_DIR
cd $ROOT_OUTPUT_DIR


# Add a line for each sample to array file
ARRAY_FILE=${ROOT_OUTPUT_DIR}/scifi_pipeline.${RUN_NAME}.${JOB_DESCR}.array_file.txt
for SAMPLE_NAME in `tail -n +2 $BARCODE_ANNOTATION | cut -d , -f 1`; do
SAMPLE_DIR=${ROOT_OUTPUT_DIR}/${SAMPLE_NAME}
echo $SAMPLE_NAME $SAMPLE_DIR >> $ARRAY_FILE
done


# Now submit job array in steps
TOTAL=`tail -n +2 $BARCODE_ANNOTATION | wc -l`
for i in `seq 0 $ARRAY_SIZE $((TOTAL - 1))`; do
ARRAY="${i}-$((i + ARRAY_SIZE - 1))"

JOB_NAME=scifi_pipeline.${RUN_NAME}.${JOB_DESCR}.${ARRAY}
JOB=${ROOT_OUTPUT_DIR}/${JOB_NAME}.sh
LOG=${ROOT_OUTPUT_DIR}/${JOB_NAME}.%a.log

echo '#!/bin/env bash' > $JOB
echo 'date' >> $JOB
echo '' >> $JOB

echo "#RUN_NAME         = ${RUN_NAME}" >> $JOB
echo "#ANNOTATION       = ${BARCODE_ANNOTATION}" >> $JOB
echo "#ROOT DIRECTORY   = ${ROOT_OUTPUT_DIR}" >> $JOB
echo "#SLURM PARAMETERS = $CPUS, $MEM, $QUEUE, $TIME, $ARRAY_SIZE" >> $JOB
echo '' >> $JOB
echo 'echo SLURM_ARRAY_TASK_ID = $SLURM_ARRAY_TASK_ID' >> $JOB
echo '' >> $JOB

# Get respective line of input
echo "ARRAY_FILE=$ARRAY_FILE" >> $JOB
echo 'readarray -t ARR < $ARRAY_FILE' >> $JOB
echo 'IFS=" " read -r -a F <<< ${ARR[$SLURM_ARRAY_TASK_ID]}' >> $JOB
echo 'SAMPLE_NAME=${F[0]}' >> $JOB
echo 'SAMPLE_DIR=${F[1]}' >> $JOB

echo '' >> $JOB

# remove tags from header
echo 'samtools view -H \' >> $JOB
echo '${SAMPLE_DIR}/${SAMPLE_NAME}.ALL.STAR.Aligned.out.bam.featureCounts.bam \' >> $JOB
echo '| grep -v "@PG" | grep -v "@RG" | grep -v "@CO" \' >> $JOB
echo '> ${SAMPLE_DIR}/${SAMPLE_NAME}.ALL.STAR.Aligned.out.bam.featureCounts.header' >> $JOB
# cat new header and body
echo '{' >> $JOB
echo '    cat ${SAMPLE_DIR}/${SAMPLE_NAME}.ALL.STAR.Aligned.out.bam.featureCounts.header &' >> $JOB
echo '    samtools view ${SAMPLE_DIR}/${SAMPLE_NAME}.ALL.STAR.Aligned.out.bam.featureCounts.bam;' >> $JOB
echo '} | \' >> $JOB
echo 'samtools sort > \' >> $JOB
echo '${SAMPLE_DIR}/${SAMPLE_NAME}.ALL.STAR.Aligned.out.bam.featureCounts.sorted.bam' >> $JOB
# index
echo 'samtools index ${SAMPLE_DIR}/${SAMPLE_NAME}.ALL.STAR.Aligned.out.bam.featureCounts.sorted.bam' >> $JOB


# remove tags from header
echo 'samtools view -H \' >> $JOB
echo '${SAMPLE_DIR}/${SAMPLE_NAME}.ALL.STAR.Aligned.out.exon.bam.featureCounts.bam \' >> $JOB
echo '| grep -v "@PG" | grep -v "@RG" | grep -v "@CO" \' >> $JOB
echo '> ${SAMPLE_DIR}/${SAMPLE_NAME}.ALL.STAR.Aligned.out.exon.bam.featureCounts.header' >> $JOB
# cat new header and body
echo '{' >> $JOB
echo '    cat ${SAMPLE_DIR}/${SAMPLE_NAME}.ALL.STAR.Aligned.out.exon.bam.featureCounts.header &' >> $JOB
echo '    samtools view ${SAMPLE_DIR}/${SAMPLE_NAME}.ALL.STAR.Aligned.out.exon.bam.featureCounts.bam;' >> $JOB
echo '} | \' >> $JOB
echo 'samtools sort > \' >> $JOB
echo '${SAMPLE_DIR}/${SAMPLE_NAME}.ALL.STAR.Aligned.out.exon.bam.featureCounts.sorted.bam' >> $JOB
# index
echo 'samtools index ${SAMPLE_DIR}/${SAMPLE_NAME}.ALL.STAR.Aligned.out.exon.bam.featureCounts.sorted.bam' >> $JOB

# remove headers
echo 'rm ${SAMPLE_DIR}/${SAMPLE_NAME}.ALL.STAR.Aligned.out.bam.featureCounts.header' >> $JOB
echo 'rm ${SAMPLE_DIR}/${SAMPLE_NAME}.ALL.STAR.Aligned.out.exon.bam.featureCounts.header' >> $JOB


echo '' >> $JOB
echo 'date' >> $JOB
echo '' >> $JOB

sbatch -J $JOB_NAME \
-o $LOG --time $TIME \
-c $CPUS --mem $MEM -p $QUEUE \
--array=$ARRAY -N 1 \
$JOB

done
