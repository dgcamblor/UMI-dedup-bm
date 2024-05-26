#!/bin/bash

#-------------------------------------------------------------------------------
# 5_variant_calling.sh
#-------------------------------------------------------------------------------
# On every deduplication result stored as a BAM file, perform variant calling
# with Mutect2, Lofreq, and VarDict (Java)
#-------------------------------------------------------------------------------

source config.conf; source utils/header.sh
shopt -s expand_aliases  # Allow aliases to be used in this script

# Activate the conda environment
eval "$(conda shell.bash hook)"; conda activate dedup_bm

if [ $PROFILE == "BRP" ]; then
    TARGET=${TARGET_BRP}
    REF_GENOME=${REF_GENOME_HG19}
    SAMPLES=(${SAMPLES_BRP[@]})
    DBSNP=${DBSNP_HG19}
    MILLS_INDEL=${MILLS_INDEL_HG19}

    data_dir="data"
    results_dir="results/brp"

elif [ $PROFILE == "UMIVAR" ]; then
    TARGET=${TARGET_UMIVAR}
    REF_GENOME=${REF_GENOME_GRCH38}
    SAMPLES=(${SAMPLES_UMIVAR[@]})
    DBSNP=${DBSNP_GRCH38}
    MILLS_INDEL=${MILLS_INDEL_GRCH38}

    data_dir="data/umivar"
    results_dir="results/umivar"

fi


header "Variant calling"

for sample in ${SAMPLES[@]}; do
    dedup_prefix=${data_dir}/${sample}/${sample}.dedup

    #---------------------------------------------------------------------------
    # Preprocessing
    #---------------------------------------------------------------------------
    if [ $PREPROCESSING -eq 1 ]; then
        mkdir -p ${data_dir}/${sample}/bqsr

        for dt in ${DEDUP_TYPES[@]}; do
            echo -e "${BLUE}Computing BQSR for ${sample} with dedup type ${dt}${NC}"

            bam_path=${dedup_prefix}_${dt}.sorted.bam
            bqsr_path=${data_dir}/${sample}/bqsr/${sample}.dedup_${dt}.sorted.bam
            bqsr_table_path=${data_dir}/${sample}/bqsr/${sample}_${dt}_recal_data.table

            echo $bqsr_table_path

            gatk BaseRecalibrator \
                -R $REF_GENOME \
                -I $bam_path \
                -O $bqsr_table_path \
                --intervals $TARGET \
                --known-sites $DBSNP \
                --known-sites $MILLS_INDEL

            echo -e "${BLUE}Applying BQSR for ${sample} with dedup type ${dt}${NC}"

            gatk ApplyBQSR \
                -R $REF_GENOME \
                -I $bam_path \
                -O $bqsr_path \
                --intervals $TARGET \
                --bqsr-recal-file $bqsr_table_path

            # Index the bam file
            samtools index $bam_path
        done
        
        echo "Changing dedup prefix to BQSR directory"
        dedup_prefix=${data_dir}/${sample}/bqsr/${sample}.dedup  # Change the dedup prefix to the BQSR directory to point to the new BAM files
    else
        echo -e "${BLUE}Skipping preprocessing for ${sample}${NC}"
    fi

    #---------------------------------------------------------------------------
    # Mutect2
    #---------------------------------------------------------------------------
    echo -e "${BBLUE}Calling variants with Mutect2 for ${sample}${NC}"

    full_results_dir=${results_dir}/${sample}/mutect2; mkdir -p $full_results_dir
    
    for dt in ${DEDUP_TYPES[@]}; do
        echo -e "${BLUE}Dedup type: ${dt}${NC}"

        bam_path=${dedup_prefix}_${dt}.sorted.bam
        vcf_path=${full_results_dir}/${sample}_${dt}_target.vcf.gz

        gatk Mutect2 \
            -R $REF_GENOME \
            -I $bam_path \
            -O $vcf_path \
            --intervals $TARGET

        # Normalize the variants
        tabix -f -p vcf $vcf_path  # Needed for bcftools norm to work
        bcftools_nfs norm -m-any --check-ref w -f $REF_GENOME -o ${full_results_dir}/${sample}_${dt}_target_norm.vcf.gz $vcf_path
        tabix -f -p vcf -f ${full_results_dir}/${sample}_${dt}_target_norm.vcf.gz
    done

    #---------------------------------------------------------------------------
    # Lofreq
    #---------------------------------------------------------------------------
    echo -e "${BBLUE}Calling variants with Lofreq for ${sample}${NC}"

    full_results_dir=${results_dir}/${sample}/lofreq; mkdir -p $full_results_dir

    for dt in ${DEDUP_TYPES[@]}; do
        echo -e "${BLUE}Dedup type: ${dt}${NC}"

        bam_path=${dedup_prefix}_${dt}.sorted.bam
        vcf_path=${full_results_dir}/${sample}_${dt}_target.vcf
        
        lofreq indelqual --dindel -f $REF_GENOME $bam_path | \
            lofreq call - \
            -f $REF_GENOME \
            -l $TARGET \
            -o $vcf_path \
            --call-indels

        # Normalize the variants
        bgzip -f $vcf_path && tabix -f -p vcf $vcf_path.gz  # Needed for bcftools norm to work
        bcftools_nfs norm -m-any --check-ref w -f $REF_GENOME -o ${full_results_dir}/${sample}_${dt}_target_norm.vcf.gz $vcf_path.gz
        tabix -f -p vcf ${full_results_dir}/${sample}_${dt}_target_norm.vcf.gz
    done

    #---------------------------------------------------------------------------
    # VarDict Java
    #---------------------------------------------------------------------------
    echo -e "${BBLUE}Calling variants with VarDict Java for ${sample}${NC}"

    full_results_dir=${results_dir}/${sample}/vardict; mkdir -p $full_results_dir

    for dt in ${DEDUP_TYPES[@]}; do
        echo -e "${BLUE}Dedup type: ${dt}${NC}"

        bam_path=${dedup_prefix}_${dt}.sorted.bam
        vcf_path=${full_results_dir}/${sample}_${dt}_target.vcf

        vardict-java \
            -G $REF_GENOME \
            -f $AF_THR \
            -N ${sample}_${dt} \
            -b $bam_path \
            -th $N_CORES \
            -c 1 -S 2 -E 3 -g 4 \
            $TARGET | \
                teststrandbias.R | \
                var2vcf_valid.pl \
                    -N ${sample}_${dt} \
                    -E \
                    -f $AF_THR > $vcf_path

        # Normalize the variants
        bgzip -f $vcf_path && tabix -f -p vcf $vcf_path.gz  # Needed for bcftools norm to work
        bcftools_nfs norm -m-any --check-ref w -f $REF_GENOME -o ${full_results_dir}/${sample}_${dt}_target_norm.vcf.gz $vcf_path.gz
        tabix -f -p vcf ${full_results_dir}/${sample}_${dt}_target_norm.vcf.gz
    done
done