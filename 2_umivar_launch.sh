#!/bin/bash

#-------------------------------------------------------------------------------
# 2_umivar_launch.sh
#-------------------------------------------------------------------------------
# Script to launch the UMIvar simulation on a set of samples specified in the
# config.conf file (SAMPLES_UMIVAR). For benchmarking purposes, after testing.
#-------------------------------------------------------------------------------

source config.conf
shopt -s expand_aliases

# Activate the conda environment
eval "$(conda shell.bash hook)"; conda activate dedup_bm


for sample in ${SAMPLES_UMIVAR[@]}; do
    raw_bam=data/umivar/${sample}/${sample}_raw.bam

    # Fix seed (-s 12) to ensure reproducibility
    python3 software/UMIvar/bin/umivar.py \
        -i $raw_bam \
        -o data/umivar/${sample}/${sample}_UV.bam \
        -t FASTQ \
        -v metadata/umivar_input.csv \
        -b metadata/umivar_target.bed \
        -e 1 \
        -s $SEED
done