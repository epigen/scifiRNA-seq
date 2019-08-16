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
    *)
          # unknown option
    ;;
esac
done

echo "RUN_NAME         = ${RUN_NAME}"
echo "FLOWCELL         = ${FLOWCELL}"
echo "NUMBER OF LANES  = ${N_LANES}"
echo "ANNOTATION       = ${BARCODE_ANNOTATION}"
echo "ROOT DIRECTORY   = ${ROOT_OUTPUT_DIR}"
echo "SLURM PARAMETERS = $CPUS, $MEM, $QUEUE, $TIME"

# Don't edit from here
mkdir -p $ROOT_OUTPUT_DIR
cd $ROOT_OUTPUT_DIR
LANES=`seq 1 $N_LANES`


# # unfortunatelly, even though STAR can output to stdout 
# # and featureCounts read from stdin, one cannot pipe them as 
# # featureCounts does not support detailed BAM output with stdin
for SAMPLE_NAME in `tail -n +2 $BARCODE_ANNOTATION | cut -d , -f 1 | grep gRNA`; do
for LANE in ${LANES[@]}; do
INPUT_BAM=/scratch/users/dbarreca/private/custom_demux/scRNA/${FLOWCELL}/${FLOWCELL}_${LANE}_samples/${FLOWCELL}_${LANE}#${SAMPLE_NAME}.bam

JOB_NAME=scifi_pipeline.${SAMPLE_NAME}.${LANE}.trim
SAMPLE_DIR=${ROOT_OUTPUT_DIR}/${SAMPLE_NAME}
mkdir -p $SAMPLE_DIR #/{logs,fastqc,mapped,barcodes,expression}
JOB=${SAMPLE_DIR}/${JOB_NAME}.sh
LOG=${SAMPLE_DIR}/${JOB_NAME}.log

OUTPUT_BAM=${SAMPLE_DIR}/${SAMPLE_NAME}.${LANE}.trimmed.bam

echo '#!/bin/env bash' > $JOB

echo "date" >> $JOB

echo "" >> $JOB
echo "python3 -u \
~/trim_grnas.py \
$INPUT_BAM \
$OUTPUT_BAM" >> $JOB

echo "" >> $JOB
echo "date" >> $JOB
echo "" >> $JOB

sbatch -J $JOB_NAME \
-o $LOG --time $TIME \
-c $CPUS --mem $MEM -p $QUEUE \
$JOB
done
done
