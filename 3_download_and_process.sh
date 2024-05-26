#!/bin/bash

#-------------------------------------------------------------------------------
# 3_download_and_process.sh
#-------------------------------------------------------------------------------
# Download and process the data:
#  - For BRP, it downloads the SRA files and converts them to fastq
#  - For UMIVAR, it only creates the Umi Inside Sequence (UIS) FASTQ file
#-------------------------------------------------------------------------------

source config.conf; source utils/header.sh
shopt -s expand_aliases

# Activate the conda environment
eval "$(conda shell.bash hook)"; conda activate dedup_bm


header "Download and process data"

#-------------------------------------------------------------------------------
# UMIVAR profile
#-------------------------------------------------------------------------------
# If using the UMIVAR profile, only create the UIS FASTQ file
#-------------------------------------------------------------------------------
if [[ $PROFILE == "UMIVAR" ]]; then
    for sample in ${SAMPLES_UMIVAR[@]}; do
        # If file R1 does not exist, skip this sample
        if [ ! -f "data/umivar/${sample}/${sample}_UV_R1.fastq" ]; then
            echo "Skipping $sample: R1 file not found"
            continue
        fi

        R1=data/umivar/${sample}/${sample}_UV_R1.fastq
        R2=data/umivar/${sample}/${sample}_UV_R2.fastq

        # Rename the files to the standard for this pipeline
        mv $R1 data/umivar/${sample}/${sample}_1.fastq; R1=data/umivar/${sample}/${sample}_1.fastq
        mv $R2 data/umivar/${sample}/${sample}_2.fastq; R2=data/umivar/${sample}/${sample}_2.fastq

        # Add the UMI to the FASTQ read (requirement for some deduplication tools)
        echo "UMI to sequence $sample"  
        python3 utils/umi2seq.py --file $R1  # Produces ${R1%.fastq}_UIS.fastq

        R1_UIS=${R1%.fastq}_UIS.fastq

        # Compress the fastq files
        echo "Compressing $sample"
        pigz -f ${R1}; R1=${R1}.gz
        pigz -f ${R2}; R2=${R2}.gz
        pigz -f ${R1_UIS}; R1_UIS=${R1_UIS}.gz
    done

    exit 0
fi

#-------------------------------------------------------------------------------
# BRP profile
#-------------------------------------------------------------------------------
# Add the sratoolkit to the PATH (if not already added)
if [[ ! $PATH =~ .*sratoolkit.* ]]; then
    export PATH=$PATH:$PWD/software/sratoolkit.3.0.7-ubuntu64/bin
fi

# For each sample in BRP_AccList.txt, convert the corresponding file to fastq using fasterq-dump
for sample in ${SAMPLES_BRP[@]}; do
    # If the pipeline.done file exists, skip this accession
    if [ -f "data/${sample}/download.done" ]; then
        echo "Skipping $sample"
        continue
    fi

    # Download the SRA file
    echo "Downloading $sample"
    prefetch $sample
    
    # Convert to fastq
    echo "Converting $sample to fastq"
    fasterq-dump "${sample}/${sample}.sra" -e $N_CORES -3 -p
    rm -r ${sample}  # Remove the SRA file

    R1=${sample}_1.fastq
    R2=${sample}_2.fastq

    mkdir -p data/${sample}

    mv $R1 data/${sample}/${R1}; R1=data/${sample}/${R1}
    mv $R2 data/${sample}/${R2}; R2=data/${sample}/${R2}

    # Unmask the UMIs in the reads
    # UMIs in the BRP data are masked, so this script creates UMIs for those
    echo "Unmask $sample"
    python3 utils/unmask_umis.py --fastq1 $R1 --fastq2 $R2 --output1 ${R1}.tmp --output2 ${R2}.tmp

    mv ${R1}.tmp $R1
    mv ${R2}.tmp $R2

    # Dictionary of masked to unmasked UMIs (output of unmask_umis.py)
    mv "${sample}.pkl" data/${sample}/

    # Add the UMI to the FASTQ read (requirement for some deduplication tools)
    echo "UMI to sequence $sample"
    python3 utils/umi2seq.py --file $R1  # Produces ${R1%.fastq}_UIS.fastq

    R1_UIS=${R1%.fastq}_UIS.fastq

    # Compress the fastq files
    echo "Compressing $sample"
    pigz -f ${R1}; R1=${R1}.gz
    pigz -f ${R2}; R2=${R2}.gz
    pigz -f ${R1_UIS}; R1_UIS=${R1_UIS}.gz

    # End the pipeline
    touch data/${sample}/download.done
done