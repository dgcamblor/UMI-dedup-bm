#!/bin/bash

#-------------------------------------------------------------------------------
# 1_test_umivar.sh
#-------------------------------------------------------------------------------
# Script to test the UMIvar simulation on a set of samples specified in the
# config.conf file (SAMPLES_UMIVAR_TEST).
# The INPUT for UMIvar is taken from undeduplicated FASTQ files, that require
# mapping to the reference genome.
#-------------------------------------------------------------------------------

source config.conf
shopt -s expand_aliases

# Activate the conda environment
eval "$(conda shell.bash hook)"; conda activate dedup_bm


for sample in "${SAMPLES_UMIVAR_TEST[@]}"; do
  echo "Processing sample: $sample"

  # Extract UMI sequences (add UMI to read name) for READ1
  umi_tools extract \
    --bc-pattern=NNNNNNNNN \
    --stdin=test_umivar/data/${sample}-OP_UMI.fastq.gz \
    --read2-in=test_umivar/data/${sample}-OP_1.fastq.gz \
    --stdout=test_umivar/data/${sample}_umi_R1.fastq.gz \
    --read2-stdout

  # Extract UMI sequences (add UMI to read name) for READ1
  umi_tools extract \
    --bc-pattern=NNNNNNNNN \
    --stdin=test_umivar/data/${sample}-OP_UMI.fastq.gz \
    --read2-in=test_umivar/data/${sample}-OP_2.fastq.gz \
    --stdout=test_umivar/data/${sample}_umi_R2.fastq.gz \
    --read2-stdout

  # Remove adapter sequences
  cutadapt \
    --cores $N_CORES \
    -a $ILLUMINA_ADAPTER \
    -A $ILLUMINA_ADAPTER \
    -m 1 \
    -o test_umivar/data/${sample}_cut_R1.fastq \
    -p test_umivar/data/${sample}_cut_R2.fastq \
    test_umivar/data/${sample}_umi_R1.fastq.gz \
    test_umivar/data/${sample}_umi_R2.fastq.gz

  # Quality filtering
  prinseq-lite.pl \
    -fastq test_umivar/data/${sample}_cut_R1.fastq \
    -fastq2 test_umivar/data/${sample}_cut_R2.fastq \
    -min_qual_mean 30 \
    -out_good test_umivar/data/filtered_${sample} \
    -out_bad test_umivar/data/bad_${sample}

  pigz test_umivar/data/filtered_${sample}_1.fastq
  pigz test_umivar/data/filtered_${sample}_2.fastq

  # Align to the reference genome, convert to bam and sort in the same step
  bwa mem -t $N_CORES -p -K 150000000 -Y \
    -R "@RG\tID:AgilentXTHS\tSM:${sample}\tPL:Illumina" \
    $REF_GENOME_GRCH38 \
    test_umivar/data/filtered_${sample}_1.fastq.gz \
    test_umivar/data/filtered_${sample}_2.fastq.gz \
    | samtools view -q $MIN_QUAL_BAM -b -h - \
    | samtools sort -@ $N_CORES -o test_umivar/data/${sample}.sorted.bam -

  samtools index test_umivar/data/${sample}.sorted.bam

  # Remove intermediate files
  rm test_umivar/data/${sample}_umi_R1.fastq.gz
  rm test_umivar/data/${sample}_umi_R2.fastq.gz

  rm test_umivar/data/${sample}_cut_R1.fastq
  rm test_umivar/data/${sample}_cut_R2.fastq

  rm test_umivar/data/filtered_${sample}_1.fastq.gz
  rm test_umivar/data/filtered_${sample}_2.fastq.gz

  rm test_umivar/data/*singletons*
  rm test_umivar/data/*bad*

done

for sample in ${SAMPLES_UMIVAR_TEST[@]}; do
    raw_bam=data/test_umivar/${sample}/${sample}.bam

    python3 software/UMIvar/bin/umivar.py \
        -i $raw_bam \
        -o data/test_umivar/${sample}/${sample}_UV.bam \
        -t FASTQ \
        -v metadata/umivar_test_input.csv \
        -b metadata/umivar_target.bed \
        -s $SEED
done

for sample in ${SAMPLES_UMIVAR_TEST[@]}; do
  # Align to the reference genome, convert to bam and sort in the same step
  bwa mem -t $N_CORES -p -K 150000000 -Y \
    -R "@RG\tID:AgilentXTHS\tSM:${sample}\tPL:Illumina" \
    $REF_GENOME_GRCH38 \
    data/test_umivar/${sample}/${sample}_UV_R1.fastq \
    data/test_umivar/${sample}/${sample}_UV_R2.fastq \
    | samtools view -q MIN_QUAL_BAM -b -h - \
    | samtools sort -@ $N_CORES -o data/test_umivar/${sample}/${sample}_UV.sorted.bam -

  samtools index data/test_umivar/${sample}/${sample}_UV.sorted.bam

done

for sample in ${SAMPLES_UMIVAR_TEST[@]}; do
    # Variant calling with VarDict
    echo "Calling with VarDict for sample: $sample"
    vardict-java \
        -G $REF_GENOME_GRCH38 \
        -f 0.005 \
        -N ${sample} \
        -b data/test_umivar/${sample}/${sample}_UV.sorted.bam \
        -th $N_CORES \
        -c 1 -S 2 -E 3 -g 4 \
        metadata/umivar_target.bed | \
            teststrandbias.R | \
            var2vcf_valid.pl \
                -N ${sample} \
                -E \
                -f 0.005 > data/test_umivar/${sample}/${sample}_vardict.vcf
done