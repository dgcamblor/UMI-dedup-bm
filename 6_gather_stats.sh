#!/bin/bash

#-------------------------------------------------------------------------------
# 6_gather_stats.sh
#-------------------------------------------------------------------------------
# Obtain different stats from the benchmarking.
# Stats are obtained from UMItools deduplication files.
#-------------------------------------------------------------------------------

source config.conf; source utils/header.sh
shopt -s expand_aliases  # Allow aliases to be used in this script

# Activate the conda environment
eval "$(conda shell.bash hook)"; conda activate dedup_bm

if [ $PROFILE == "BRP" ]; then
    SAMPLES=(${SAMPLES_BRP[@]})

    TARGET="metadata/LBx_BRP_hg19.bed"
    DATA_DIR="data"
    STATS_DIR="stats"

elif [ $PROFILE == "UMIVAR" ]; then
    SAMPLES=(${SAMPLES_UMIVAR[@]})

    TARGET="metadata/umivar_target.bed"
    DATA_DIR="data/umivar"
    STATS_DIR="stats/umivar"

fi


header "Gathering stats"

echo "sample,umi_ratio,median_cov" > ${STATS_DIR}/stats_summary.csv

for sample in ${SAMPLES[@]}; do
    echo "Processing sample: $sample"
    # umi_ratio-----------------------------------------------------------------
    umi_ratio=$(awk 'NR>1 {sum += $1 * $3; total += $3} END {print sum / total}' ${STATS_DIR}/${sample}/${sample}.dedup_UT.stats_per_umi_per_position.tsv)

    # median_cov----------------------------------------------------------------
    echo "mosdepth for sample: $sample"

    if [ $PROFILE == "BRP" ]; then
        # Create on target
        echo "Creating on-target for per-base cov computing"
        samtools view -b -L ${TARGET} ${DATA_DIR}/${sample}/${sample}.dedup_UT.sorted.bam > ${DATA_DIR}/${sample}/${sample}.dedup_UT.sorted_on_target.bam
        samtools index ${DATA_DIR}/${sample}/${sample}.dedup_UT.sorted_on_target.bam
        input_bam=${DATA_DIR}/${sample}/${sample}.dedup_UT.sorted_on_target.bam

    elif [ $PROFILE == "UMIVAR" ]; then
        input_bam=${DATA_DIR}/${sample}/${sample}.dedup_UT.sorted.bam
    
    fi

    echo -e "${BBLUE}Computing per-base coverage for ${sample}${NC}"
    mosdepth \
        --threads 30 \
        --use-median \
        --by ${TARGET} \
        ${sample} \
        ${input_bam}

    if [ $PROFILE == "BRP" ]; then
        rm ${DATA_DIR}/${sample}/${sample}.dedup_UT.sorted_on_target.bam
        rm ${DATA_DIR}/${sample}/${sample}.dedup_UT.sorted_on_target.bam.bai
    fi

    mv ${sample}* ${STATS_DIR}/${sample}

    median_cov=$(zcat ${STATS_DIR}/${sample}/${sample}.per-base.bed.gz | datamash median 4)

    echo "${sample},${umi_ratio},${median_cov}" >> ${STATS_DIR}/stats_summary.csv
done
