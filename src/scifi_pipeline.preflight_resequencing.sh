
# The goal is to make the structure of demultiplexed files for one experiment
# that was resequenced look like it is one sequencing run, simply with
# additional lanes.
# Since the original files are in a location non-writable for me, the new
# structure will be created in a new location too.

# For PD190:
# sh src/scifi_pipeline.preflight_resequencing.sh \
#   BSF_0666_HCW2YDRXX 2 \
#   BSF_0670_HFCK2DRXX 2
# RUN_1=BSF_0666_HCW2YDRXX
# RUN_1_LANES=2
# RUN_2=BSF_0670_HFCK2DRXX
# RUN_2_LANES=2

# For PD193_fivelines_383k and PD193_humanmouse_765k:
# sh src/scifi_pipeline.preflight_resequencing.sh \
#   BSF_0674_HL7L5DMXX 2 \
#   BSF_0678_HM7GMDMXX 2
# RUN_1=BSF_0674_HL7L5DMXX
# RUN_1_LANES=2
# RUN_2=BSF_0678_HM7GMDMXX
# RUN_2_LANES=2

# For PD195-1_Tcells_765k-sample1_P7:
# sh src/scifi_pipeline.preflight_resequencing.sh \
#   BSF_0684_HGNMYDRXX 2 \
#   BSF_0688_HMCGTDMXX 2
# RUN_1=BSF_0684_HGNMYDRXX
# RUN_1_LANES=2
# RUN_2=BSF_0688_HMCGTDMXX
# RUN_2_LANES=2

RUN_1=$1
RUN_1_LANES=$2

RUN_2=$3
RUN_2_LANES=$4

ORIGINAL_ROOT_DIR=/scratch/users/dbarreca/private/custom_demux/scRNA/
NEW_ROOT_DIR=/scratch/users/${USER}

mkdir -p $NEW_ROOT_DIR/custom_demux/scRNA/
cd $NEW_ROOT_DIR/custom_demux/scRNA/


NEW_RUN_NAME=${RUN_1}
NEW_TOTAL_LANES=$((RUN_1_LANES + RUN_2_LANES))

I=1
for N in `seq $RUN_1_LANES`; do
mkdir -p ${NEW_RUN_NAME}_${I}_samples

cd ${NEW_RUN_NAME}_${I}_samples
ln -s ${ORIGINAL_ROOT_DIR}/${RUN_1}/${RUN_1}_${I}_samples/${RUN_1}_${I}#*.bam .
cd ..
I=$((I + 1))
done

for N in `seq $RUN_2_LANES`; do
mkdir -p ${NEW_RUN_NAME}_${I}_samples

cd ${NEW_RUN_NAME}_${I}_samples
ln -s ${ORIGINAL_ROOT_DIR}/${RUN_2}/${RUN_2}_${N}_samples/${RUN_2}_${N}#*.bam .
rename ${RUN_2}_${N}# ${RUN_1}_${I}# *.bam
cd ..

I=$((I + 1))
done


if [ $((I - 1)) == $NEW_TOTAL_LANES ]; then
    echo 'SUCCESS!'
else
    echo 'ERROR!'
fi
