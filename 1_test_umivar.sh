#!/bin/bash

#-------------------------------------------------------------------------------
# 1_test_umivar.sh
#-------------------------------------------------------------------------------
# Script to test the UMIvar simulation on a set of samples specified in the
# config.conf file (SAMPLES_UMIVAR_TEST).
# The input for UMIvar is taken from undeduplicated FASTQ files, which require
# processing and mapping to the reference genome.
#-------------------------------------------------------------------------------

source config.conf
shopt -s expand_aliases

# Activate the conda environment
eval "$(conda shell.bash hook)"; conda activate dedup_bm


data_dir="data/umivar_test"
results_dir="results/umivar_test"

for sample in "${SAMPLES_UMIVAR_TEST[@]}"; do
  echo "Processing sample: $sample"

  # Extract UMI sequences (add UMI to read name) for READ1
  umi_tools extract \
    --bc-pattern=NNNNNNNNN \
    --stdin=${data_dir}/${sample}/${sample}_UMI.fastq.gz \
    --read2-in=${data_dir}/${sample}/${sample}_1.fastq.gz \
    --stdout=${data_dir}/${sample}/${sample}_umi_R1.fastq.gz \
    --read2-stdout

  # Extract UMI sequences (add UMI to read name) for READ1
  umi_tools extract \
    --bc-pattern=NNNNNNNNN \
    --stdin=${data_dir}/${sample}/${sample}_UMI.fastq.gz \
    --read2-in=${data_dir}/${sample}/${sample}_2.fastq.gz \
    --stdout=${data_dir}/${sample}/${sample}_umi_R2.fastq.gz \
    --read2-stdout

  # Remove adapter sequences
  cutadapt \
    --cores $N_CORES \
    -a $ILLUMINA_ADAPTER \
    -A $ILLUMINA_ADAPTER \
    -m 1 \
    -o ${data_dir}/${sample}/${sample}_cut_R1.fastq \
    -p ${data_dir}/${sample}/${sample}_cut_R2.fastq \
    ${data_dir}/${sample}/${sample}_umi_R1.fastq.gz \
    ${data_dir}/${sample}/${sample}_umi_R2.fastq.gz

  # Quality filtering
  prinseq-lite.pl \
    -fastq ${data_dir}/${sample}/${sample}_cut_R1.fastq \
    -fastq2 ${data_dir}/${sample}/${sample}_cut_R2.fastq \
    -min_qual_mean 30 \
    -out_good ${data_dir}/${sample}/filtered_${sample} \
    -out_bad ${data_dir}/${sample}/bad_${sample}

  pigz ${data_dir}/${sample}/filtered_${sample}_1.fastq
  pigz ${data_dir}/${sample}/filtered_${sample}_2.fastq

done

for sample in "${SAMPLES_UMIVAR_TEST[@]}"; do
  echo "Aligning sample: $sample"
  ${bwa2} mem -t $N_CORES -K 150000000 -Y \
    -R "@RG\tID:AgilentXTHS\tSM:${sample}\tPL:Illumina" \
    $REF_GENOME_GRCH38 \
    ${data_dir}/${sample}/filtered_${sample}_1.fastq.gz \
    ${data_dir}/${sample}/filtered_${sample}_2.fastq.gz \
    | samtools view -q $MIN_QUAL_BAM -b -h - \
    | samtools sort -@ $N_CORES -o ${data_dir}/${sample}/${sample}.sorted.bam -

  samtools index ${data_dir}/${sample}/${sample}.sorted.bam

  # Remove intermediate files
  rm ${data_dir}/${sample}/${sample}_umi_R1.fastq.gz
  rm ${data_dir}/${sample}/${sample}_umi_R2.fastq.gz

  rm ${data_dir}/${sample}/${sample}_cut_R1.fastq
  rm ${data_dir}/${sample}/${sample}_cut_R2.fastq

  rm ${data_dir}/${sample}/filtered_${sample}_1.fastq.gz
  rm ${data_dir}/${sample}/filtered_${sample}_2.fastq.gz

  rm ${data_dir}/${sample}/*singletons*
  rm ${data_dir}/${sample}/*bad*
done

# Simulate variants with UMIvar
for sample in ${SAMPLES_UMIVAR_TEST[@]}; do
    raw_bam=${data_dir}/${sample}/${sample}.sorted.bam

    python3 software/UMIvar/bin/umivar.py \
        -i $raw_bam \
        -o ${data_dir}/${sample}/${sample}_UV.bam \
        -t FASTQ \
        -v metadata/umivar_test_input.csv \
        -b metadata/umivar_target.bed \
        -e 1 \
        -s $SEED
done

# Align the UMIvar output to the reference genome
for sample in ${SAMPLES_UMIVAR_TEST[@]}; do
  # Align to the reference genome, convert to bam and sort in the same step
  ${bwa2} mem -t $N_CORES -K 150000000 -Y \
    -R "@RG\tID:AgilentXTHS\tSM:${sample}\tPL:Illumina" \
    $REF_GENOME_GRCH38 \
    ${data_dir}/${sample}/${sample}_UV_R1.fastq \
    ${data_dir}/${sample}/${sample}_UV_R2.fastq \
    | samtools view -q ${MIN_QUAL_BAM} -b -h - \
    | samtools sort -@ ${N_CORES} -o ${data_dir}/${sample}/${sample}_UV.sorted.bam -

  samtools index ${data_dir}/${sample}/${sample}_UV.sorted.bam

done

# Deduplicate the UMIvar output with UMI-tools
for sample in ${SAMPLES_UMIVAR_TEST[@]}; do
  umi_tools dedup \
    -I ${data_dir}/${sample}/${sample}_UV.sorted.bam \
    -S ${data_dir}/${sample}/${sample}_UV.dedup.sorted.bam \
    --paired

  samtools index ${data_dir}/${sample}/${sample}_UV.dedup.sorted.bam
done

mkdir -p ${results_dir}

for sample in ${SAMPLES_UMIVAR_TEST[@]}; do
    # Variant calling with VarDict
    echo "Calling with VarDict for sample: $sample"
    vardict-java \
        -G $REF_GENOME_GRCH38 \
        -f 0.005 \
        -N ${sample} \
        -b ${data_dir}/${sample}/${sample}_UV.dedup.sorted.bam \
        -th $N_CORES \
        -c 1 -S 2 -E 3 -g 4 \
        metadata/umivar_target.bed | \
            teststrandbias.R | \
            var2vcf_valid.pl \
                -N ${sample} \
                -E \
                -f 0.005 > ${results_dir}/${sample}_vardict.vcf

    # Normalize the variants
    echo -e "${BBLUE}Normalizing the variants for ${sample}${NC}"
    bgzip -f ${results_dir}/${sample}_vardict.vcf && tabix -f -p vcf ${results_dir}/${sample}_vardict.vcf.gz
    bcftools norm -m-any --check-ref w -f $REF_GENOME_GRCH38 -o ${results_dir}/${sample}_vardict_norm.vcf.gz ${results_dir}/${sample}_vardict.vcf.gz
    
    # Unzip the normalized VCF (required for comp_umivar_test.py)
    gunzip ${results_dir}/${sample}_vardict_norm.vcf.gz
done