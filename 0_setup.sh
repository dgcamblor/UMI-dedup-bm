#!/bin/bash

#-------------------------------------------------------------------------------
# 0_setup.sh
#-------------------------------------------------------------------------------
# Configure the setup (software) to perform the benchmarking on deduplication
# software.
#-------------------------------------------------------------------------------

mkdir -p software data 

source config.conf

# Create the conda environment
conda create -n dedup_bm python=3.10 -y
eval "$(conda shell.bash hook)"; conda activate dedup_bm


#-------------------------------------------------------------------------------
# UMIvar
#-------------------------------------------------------------------------------

conda install bioconda::pysam -y
git clone https://github.com/dgcamblor/UMIvar
mv UMIvar software/

#-------------------------------------------------------------------------------
# General tools installation
#-------------------------------------------------------------------------------

# Picard
wget https://github.com/broadinstitute/picard/releases/download/1.131/picard-tools-1.131.zip
mv picard-tools-1.131.zip software/
unzip software/picard-tools-1.131.zip -d software/
rm software/picard-tools-1.131.zip

# Datamash
conda install datamash -y

# Dependencies
conda install -c conda-forge openjdk samtools ncurses -y

# SRA toolkit
wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -vxzf sratoolkit.tar.gz
mv sratoolkit.* software/
rm software/sratoolkit.tar.gz

# BWA
conda install -c bioconda bwa -y

# bedtools
conda install -c bioconda bedtools -y

# bcftools
conda install -c bioconda bcftools -y

# tabix
conda install -c bioconda tabix -y

#-------------------------------------------------------------------------------
# Variant calling
#-------------------------------------------------------------------------------

# GATK / Picard (+ Mutect2)
conda install -c bioconda gatk4 -y

# Lofreq
conda install -c bioconda lofreq -y

# VarDict
conda install -c bioconda vardict-java -y

#-------------------------------------------------------------------------------
# Deduplication tools
#-------------------------------------------------------------------------------

# UMI-tools
pip install umi_tools==1.1.4  # conda install -c bioconda umi_tools=1.1.4

# fgbio
wget https://github.com/fulcrumgenomics/fgbio/releases/download/2.1.0/fgbio-2.1.0.jar
mkdir -p software/fgbio
mv fgbio-2.1.0.jar software/fgbio/

# UMICollapse
conda install -c bioconda umicollapse -y

#! UMIc  #! Not considered because it does not work with very deep sequencing data
# conda install -c r r-essentials -y
# conda install -c conda-forge r-tidyverse r-data.table r-stringdist r-pryr -y  # Dependencies
# conda install -c bioconda bioconductor-Biostrings bioconductor-ShortRead -y  # Dependencies
# git clone https://github.com/BiodataAnalysisGroup/UMIc
# mv UMIc/ software/

# gencore
conda install -c bioconda gencore -y

#-------------------------------------------------------------------------------
# File processing
#-------------------------------------------------------------------------------

# Create on-target true variants (SEQC2 known positives and BRP sequenced regions)
bedtools intersect -a metadata/KnownPositives_hg19.vcf.gz -b metadata/LBx_BRP_hg19.bed -header > metadata/KnownPositives_hg19_target.vcf

# Normalize VCF
bcftools norm -m-any -o metadata/KnownPositives_hg19_target_norm.vcf metadata/KnownPositives_hg19_target.vcf