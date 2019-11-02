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
    --expected-cell-number=*)
    EXPECTED_CELL_NUMBER="${i#*=}"
    shift # past argument=value
    ;;
    --min-umi-output=*)
    MIN_UMI_OUTPUT="${i#*=}"
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
    --correct-barcodes=*)
    CORRECT_BARCODES="${i#*=}"
    shift # past argument=value
    ;;
    --correct-barcode-file=*)
    CORRECT_BARCODE_FILE="${i#*=}"
    shift # past argument=value
    ;;
    --no-overwrite=*)
    NO_OVERWRITE="${i#*=}"
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
echo "VARIABLES            = ${VARIABLES}"
echo "SPECIES_MIXTURE      = ${SPECIES_MIXTURE}"
echo "ROOT DIRECTORY       = ${ROOT_OUTPUT_DIR}"
echo "SLURM PARAMETERS     = $CPUS, $MEM, $QUEUE, $TIME, $ARRAY_SIZE"
echo "CORRECT_BARCODES     = $CORRECT_BARCODES"
echo "CORRECT_BARCODE_FILE = $CORRECT_BARCODE_FILE"

# Start

# # Make path absolute
if [[ ! "$BARCODE_ANNOTATION" = /* ]]; then
    BARCODE_ANNOTATION=`pwd`/${BARCODE_ANNOTATION}
fi

SUMMARIZER=`pwd`/src/scifi_pipeline.summarizer.py


ADDITIONAL_ARGS=""
JOB_DESCR="filter"
if [[ $SPECIES_MIXTURE = "1" ]]; then
    ADDITIONAL_ARGS="$ADDITIONAL_ARGS --species-mixture "
fi
if [[ $CORRECT_BARCODES = "1" ]]; then
    ADDITIONAL_ARGS="$ADDITIONAL_ARGS --correct-r2-barcodes "
    JOB_DESCR="filter_corrected"
fi
if [[ ! -z $CORRECT_BARCODE_FILE ]]; then
    ADDITIONAL_ARGS="$ADDITIONAL_ARGS --correct-r2-barcode-file $CORRECT_BARCODE_FILE"
fi

if [[ ! -z NO_OVERWRITE ]]; then
    NO_OVERWRITE="1"
fi

mkdir -p $ROOT_OUTPUT_DIR
cd $ROOT_OUTPUT_DIR


# Add a line for each sample to array file
ARRAY_FILE=${ROOT_OUTPUT_DIR}/scifi_pipeline.${RUN_NAME}.${JOB_DESCR}.array_file.txt
rm -rf $ARRAY_FILE
for SAMPLE_NAME in `tail -n +2 $BARCODE_ANNOTATION | cut -d , -f 1`; do
    SAMPLE_DIR=${ROOT_OUTPUT_DIR}/${SAMPLE_NAME}
    if [[ $NO_OVERWRITE = "1" ]]; then
        if [[ ! -f ${SAMPLE_DIR}/${SAMPLE_NAME}.metrics.csv.gz ]]; then
            echo $SAMPLE_NAME $SAMPLE_DIR >> $ARRAY_FILE
        fi
    else
        echo $SAMPLE_NAME $SAMPLE_DIR >> $ARRAY_FILE
    fi
done


# Now submit job array in steps
TOTAL=`cat $ARRAY_FILE | wc -l`

# # reduce array size if not enough samples to do
TOTAL=$((TOTAL<ARRAY_SIZE ? TOTAL : ARRAY_SIZE))
ARRAY_SIZE=$((TOTAL<ARRAY_SIZE ? TOTAL : ARRAY_SIZE))

for i in `seq 0 $ARRAY_SIZE $((TOTAL - 1))`; do
ARRAY="${i}-$((i + ARRAY_SIZE - 1))"

JOB_NAME=scifi_pipeline.${RUN_NAME}.${JOB_DESCR}.${ARRAY}
JOB=${ROOT_OUTPUT_DIR}/${JOB_NAME}.sh
LOG=${ROOT_OUTPUT_DIR}/${JOB_NAME}.%a.log

echo '#!/bin/env bash' > $JOB
echo 'date' >> $JOB
echo '' >> $JOB

echo "#RUN_NAME         = ${RUN_NAME}" >> $JOB
echo "#FLOWCELL         = ${FLOWCELL}" >> $JOB
echo "#NUMBER OF LANES  = ${N_LANES}" >> $JOB
echo "#BARCODE NUMBER   = ${N_BARCODES}" >> $JOB
echo "#ANNOTATION       = ${BARCODE_ANNOTATION}" >> $JOB
echo "#ROOT DIRECTORY   = ${ROOT_OUTPUT_DIR}" >> $JOB
echo "#STAR EXECUTABLE  = ${STAR_EXE}" >> $JOB
echo "#STAR DIRECTORY   = ${STAR_DIR}" >> $JOB
echo "#GTF FILE         = ${GTF_FILE}" >> $JOB
echo "#SLURM PARAMETERS = $CPUS, $MEM, $QUEUE, $TIME, $ARRAY_SIZE" >> $JOB
echo '' >> $JOB
echo 'echo SLURM_ARRAY_TASK_ID = $SLURM_ARRAY_TASK_ID' >> $JOB
echo '' >> $JOB

# Get respective line of input
echo "ARRAY_FILE=$ARRAY_FILE" >> $JOB
echo 'readarray -t ARR < $ARRAY_FILE' >> $JOB
echo 'IFS=" " read -r -a F <<< ${ARR[$SLURM_ARRAY_TASK_ID]}' >> $JOB
echo 'SAMPLE_NAME=${F[0]}' >> $JOB
echo 'echo $SAMPLE_NAME' >> $JOB
echo 'SAMPLE_DIR=${F[1]}' >> $JOB
echo 'echo $SAMPLE_DIR' >> $JOB

echo '' >> $JOB
echo "python3 -u $SUMMARIZER \\
--r1-annot $BARCODE_ANNOTATION \\
--r1-attributes $VARIABLES \\
--cell-barcodes r2 \\
--only-summary \\
--no-save-intermediate \\
--min-umi-output $MIN_UMI_OUTPUT \\
--expected-cell-number $EXPECTED_CELL_NUMBER \\
--no-output-header \\
--save-gene-expression \\
--correct-r1-barcodes \\
$ADDITIONAL_ARGS \\" >> $JOB
echo '--sample-name $SAMPLE_NAME \
${SAMPLE_DIR}/${SAMPLE_NAME}.*.STAR.Aligned.out.bam.featureCounts.bam \
${SAMPLE_DIR}/${SAMPLE_NAME}' >> $JOB
echo '' >> $JOB

echo 'date' >> $JOB
echo '' >> $JOB
echo "python3 -u $SUMMARIZER \\
--r1-annot $BARCODE_ANNOTATION \\
--r1-attributes $VARIABLES \\
--cell-barcodes r2 \\
--only-summary \\
--no-save-intermediate \\
--min-umi-output $MIN_UMI_OUTPUT \\
--expected-cell-number $EXPECTED_CELL_NUMBER \\
--no-output-header \\
--save-gene-expression \\
--correct-r1-barcodes \\
$ADDITIONAL_ARGS \\" >> $JOB
echo '--sample-name $SAMPLE_NAME \
${SAMPLE_DIR}/${SAMPLE_NAME}.*.STAR.Aligned.out.exon.bam.featureCounts.bam \
${SAMPLE_DIR}/${SAMPLE_NAME}.exon' >> $JOB
echo '' >> $JOB
echo 'date' >> $JOB
echo '' >> $JOB

sbatch -J $JOB_NAME \
-o $LOG --time $TIME \
-c $CPUS --mem $MEM -p $QUEUE \
--array=$ARRAY -N 1 \
$JOB

done
