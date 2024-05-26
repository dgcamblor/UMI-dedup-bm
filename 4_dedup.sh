#!/bin/bash

#-------------------------------------------------------------------------------
# 4_dedup.sh
#-------------------------------------------------------------------------------
# Perform deduplication using different tools
#-------------------------------------------------------------------------------

source config.conf; source utils/header.sh
shopt -s expand_aliases  # Allow aliases to be used in this script

# Activate the conda environment
eval "$(conda shell.bash hook)"; conda activate dedup_bm

if [ $PROFILE == "BRP" ]; then
    TARGET=${TARGET_BRP}
    REF_GENOME=${REF_GENOME_HG19}
    SAMPLES=(${SAMPLES_BRP[@]})
    READ_STRUCTURE=${READ_STRUCTURE_BRP}
    MAX_EDITS=0  # Force identity (UMI sequences are arbitrarily assigned to masked UMIs)

    data_dir="data"
    stats_dir="stats"

elif [ $PROFILE == "UMIVAR" ]; then
    TARGET=${TARGET_UMIVAR}
    REF_GENOME=${REF_GENOME_GRCH38}
    SAMPLES=(${SAMPLES_UMIVAR[@]})
    READ_STRUCTURE=${READ_STRUCTURE_UMIVAR}
    MAX_EDITS=1  # Allow 1 edit distance for UMI sequences (usual default)

    data_dir="data/umivar"
    stats_dir="stats/umivar"

fi


header "Deduplication"

for sample in ${SAMPLES[@]}; do
    R1=${data_dir}/${sample}/${sample}_1.fastq.gz
    R2=${data_dir}/${sample}/${sample}_2.fastq.gz
    R1_UIS=${data_dir}/${sample}/${sample}_1_UIS.fastq.gz

    mkdir -p ${stats_dir}/${sample}

    # Align the reads
    echo -e "${BBLUE}Aligning reads for ${sample}, creating sorted BAM${NC}"
    ${bwa2} mem -t $N_CORES -K 150000000 -Y -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" $REF_GENOME $R1 $R2 | \
        samtools view -b - | \
        samtools sort -@ $N_CORES -o ${data_dir}/${sample}/${sample}.sorted.bam -

    if [ $PREPROCESSING -eq 1 ]; then
        echo -e "${BLUE}Filtering reads with mapping quality < $MIN_MAPQ for ${sample}${NC}"  # Default: 20
        samtools view -b -q $MIN_MAPQ ${data_dir}/${sample}/${sample}.sorted.bam > ${data_dir}/${sample}/${sample}.sorted.filtered.bam

        mv ${data_dir}/${sample}/${sample}.sorted.bam ${data_dir}/${sample}/${sample}.sorted.unfiltered.bam
        mv ${data_dir}/${sample}/${sample}.sorted.filtered.bam ${data_dir}/${sample}/${sample}.sorted.bam
    fi

    # Sort and index
    echo -e "${BBLUE}Indexing BAM for ${sample}${NC}"
    samtools index ${data_dir}/${sample}/${sample}.sorted.bam  # UMI in header, sep: _

    #---------------------------------------------------------------------------
    # GATK Picard MarkDuplicates
    #---------------------------------------------------------------------------

    # Without BARCODE_TAG ------------------------------------------------------

    echo -e "${BBLUE}Marking duplicates for ${sample} (BARCODE_TAG disabled)${NC}"
    gatk MarkDuplicates \
        -I ${data_dir}/${sample}/${sample}.sorted.bam \
        -O ${data_dir}/${sample}/${sample}.dedup_MD.bam \
        -M ${stats_dir}/${sample}/${sample}.dedup_MD.metrics \
        --REMOVE_DUPLICATES true 

    # Sort and index
    samtools sort -@ $N_CORES -o ${data_dir}/${sample}/${sample}.dedup_MD.sorted.bam ${data_dir}/${sample}/${sample}.dedup_MD.bam; rm ${data_dir}/${sample}/${sample}.dedup_MD.bam
    samtools index ${data_dir}/${sample}/${sample}.dedup_MD.sorted.bam

    # With BARCODE_TAG ---------------------------------------------------------

    # Start by adding the UMI to the BAM file (RX tag)
    python3 utils/umi2tag.py --input ${data_dir}/${sample}/${sample}.sorted.bam --output ${data_dir}/${sample}/${sample}.sorted.tags.bam

    gatk MarkDuplicates \
        -I ${data_dir}/${sample}/${sample}.sorted.tags.bam \
        -O ${data_dir}/${sample}/${sample}.dedup_MDBT.bam \
        -M ${stats_dir}/${sample}/${sample}.dedup_MDBT.metrics \
        --BARCODE_TAG RX \
        --REMOVE_DUPLICATES true 

    # Sort and index
    samtools sort -@ $N_CORES -o ${data_dir}/${sample}/${sample}.dedup_MDBT.sorted.bam ${data_dir}/${sample}/${sample}.dedup_MDBT.bam; rm ${data_dir}/${sample}/${sample}.dedup_MDBT.bam
    samtools index ${data_dir}/${sample}/${sample}.dedup_MDBT.sorted.bam

    #---------------------------------------------------------------------------
    # UMI-tools dedup
    #---------------------------------------------------------------------------

    # If using BRP profile, use UMI identities without any edits
    if [ $PROFILE == "BRP" ]; then
        echo -e "${BBLUE}Deduplicating with UMI-tools for ${sample} (no edits)${NC}"
        umi_tools dedup \
            -I ${data_dir}/${sample}/${sample}.sorted.bam \
            --paired \
            -S ${data_dir}/${sample}/${sample}.dedup_UT.bam \
            --method unique \
            --output-stats ${stats_dir}/${sample}/${sample}.dedup_UT.stats \
            --no-sort-output

    elif [ $PROFILE == "UMIVAR" ]; then
        echo -e "${BBLUE}Deduplicating with UMI-tools for ${sample} (edits = 1)${NC}"
        umi_tools dedup \
            -I ${data_dir}/${sample}/${sample}.sorted.bam \
            --paired \
            -S ${data_dir}/${sample}/${sample}.dedup_UT.bam \
            --method adjacency \
            --edit-distance-threshold $MAX_EDITS \
            --output-stats ${stats_dir}/${sample}/${sample}.dedup_UT.stats \
            --no-sort-output

    fi
    
    # Sort and index
    samtools sort -@ $N_CORES -o ${data_dir}/${sample}/${sample}.dedup_UT.sorted.bam ${data_dir}/${sample}/${sample}.dedup_UT.bam; rm ${data_dir}/${sample}/${sample}.dedup_UT.bam
    samtools index ${data_dir}/${sample}/${sample}.dedup_UT.sorted.bam

    #---------------------------------------------------------------------------
    # fgbio (best practices implementation)
    #---------------------------------------------------------------------------

    echo -e "${BBLUE}Deduplicating with fgbio for ${sample}${NC}"

    # Process the reads into an unmapped BAM file: UMIs are added to tags
    echo -e "${BLUE}Creating unmapped bam with UMIs ${sample}${NC}"
    fgbio --tmp-dir ${TMP_DIR} --compression 1 --async-io FastqToBam \
        --input ${R1_UIS} ${R2} \
        --read-structures ${READ_STRUCTURE} \
        --sample ${sample} \
        --library ${sample} \
        --platform-unit ${sample} \
        --output ${data_dir}/${sample}/${sample}.unmapped.bam

    # Convert the unmapped BAM file to a mapped BAM file and realign it, carrying the information of the unmapped BAM file
    echo -e "${BLUE}Realigning the unmapped BAM for ${sample}${NC}"
    samtools fastq ${data_dir}/${sample}/${sample}.unmapped.bam \
        | ${bwa2} mem -t ${N_CORES} -p -K 150000000 -Y ${REF_GENOME} - \
        | fgbio --tmp-dir ${TMP_DIR} --compression 1 --async-io ZipperBams \
            --unmapped ${data_dir}/${sample}/${sample}.unmapped.bam \
            --ref ${REF_GENOME} \
            --output ${data_dir}/${sample}/${sample}.mapped.bam

    rm ${data_dir}/${sample}/${sample}.unmapped.bam

    # Group reads by UMI
    if [ $PROFILE == "BRP" ]; then
        echo -e "${BLUE}Grouping reads by UMI for ${sample}${NC}"
        fgbio --tmp-dir ${TMP_DIR} --compression 1 --async-io GroupReadsByUmi \
            --input ${data_dir}/${sample}/${sample}.mapped.bam \
            --strategy Identity \
            --edits $MAX_EDITS \
            --output ${data_dir}/${sample}/${sample}.grouped.bam \
            --family-size-histogram ${stats_dir}/${sample}/${sample}.dedup_FGB.grouped.hist.txt

    elif [ $PROFILE == "UMIVAR" ]; then
        fgbio --tmp-dir ${TMP_DIR} --compression 1 --async-io GroupReadsByUmi \
            --input ${data_dir}/${sample}/${sample}.mapped.bam \
            --strategy adjacency \
            --edits $MAX_EDITS \
            --output ${data_dir}/${sample}/${sample}.grouped.bam \
            --family-size-histogram ${stats_dir}/${sample}/${sample}.dedup_FGB.grouped.hist.txt

    fi

    rm ${data_dir}/${sample}/${sample}.mapped.bam

    # Call molecular consensus reads
    echo -e "${BLUE}Calling molecular consensus reads for ${sample}${NC}"
    fgbio --tmp-dir ${TMP_DIR} --compression 1 CallMolecularConsensusReads \
        --input ${data_dir}/${sample}/${sample}.grouped.bam \
        --output ${data_dir}/${sample}/${sample}.cons.unmapped.bam \
        --min-reads $SUPP_READS \
        --threads ${N_CORES}

    rm ${data_dir}/${sample}/${sample}.grouped.bam

    # Re-align the consensus reads
    echo -e "${BLUE}Re-aligning the consensus reads for ${sample}${NC}"
    samtools fastq ${data_dir}/${sample}/${sample}.cons.unmapped.bam \
        | ${bwa2} mem -t ${N_CORES} -p -K 150000000 -Y ${REF_GENOME} - \
        | fgbio --tmp-dir ${TMP_DIR} --compression 1 --async-io ZipperBams \
            --unmapped ${data_dir}/${sample}/${sample}.cons.unmapped.bam \
            --ref ${REF_GENOME} \
            --tags-to-reverse Consensus \
            --tags-to-revcomp Consensus \
            --output ${data_dir}/${sample}/${sample}.dedup_FGB.bam

    rm ${data_dir}/${sample}/${sample}.cons.unmapped.bam

    # Sort and index
    samtools sort -@ $N_CORES -o ${data_dir}/${sample}/${sample}.dedup_FGB.sorted.bam ${data_dir}/${sample}/${sample}.dedup_FGB.bam; rm ${data_dir}/${sample}/${sample}.dedup_FGB.bam
    samtools index ${data_dir}/${sample}/${sample}.dedup_FGB.sorted.bam

    #---------------------------------------------------------------------------
    # UMICollapse
    #---------------------------------------------------------------------------

    echo -e "${BBLUE}Deduplicating with UMIcollapse for ${sample}${NC}"

    umicollapse bam \
        -i ${data_dir}/${sample}/${sample}.sorted.bam \
        -o ${data_dir}/${sample}/${sample}.dedup_UCOL.bam \
        --umi-sep _ \
        --paired \
        --algo adj \
        -k $MAX_EDITS \
        --two-pass

    # Sort and index
    samtools sort -@ $N_CORES -o ${data_dir}/${sample}/${sample}.dedup_UCOL.sorted.bam ${data_dir}/${sample}/${sample}.dedup_UCOL.bam; rm ${data_dir}/${sample}/${sample}.dedup_UCOL.bam
    samtools index ${data_dir}/${sample}/${sample}.dedup_UCOL.sorted.bam

    #---------------------------------------------------------------------------
    # UMIc  #! Currently skipped due to very slow performance in BRP
    #---------------------------------------------------------------------------

    #! mkdir -p ${data_dir}/${sample}/UMIc
    #! cp ${R1_UIS} ${data_dir}/${sample}/UMIc/${sample}_R1.fastq.gz
    #! cp ${R2} ${data_dir}/${sample}/UMIc/${sample}_R2.fastq.gz

    #---------------------------------------------------------------------------
    # gencore
    #---------------------------------------------------------------------------

    # Convert all "_" to ":" in the BAM file (gencore requirement)
    samtools view -h ${data_dir}/${sample}/${sample}.sorted.bam | sed 's/_/:/g' | samtools view -b - > ${data_dir}/${sample}/${sample}.sorted.gc.bam

    echo -e "${BBLUE}Deduplicating with gencore for ${sample}${NC}"
    gencore \
        -i ${data_dir}/${sample}/${sample}.sorted.gc.bam \
        -o ${data_dir}/${sample}/${sample}.dedup_GC.bam \
        -r $REF_GENOME \
        --umi_prefix "" \
        --no_duplex \
        --supporting_reads $SUPP_READS \
        --umi_diff_threshold $MAX_EDITS \
        --json ${stats_dir}/${sample}/${sample}.dedup_GC.json \
        --html ${stats_dir}/${sample}/${sample}.dedup_GC.html

    rm ${data_dir}/${sample}/${sample}.sorted.gc.bam

    # Sort and index
    samtools sort -@ $N_CORES -o ${data_dir}/${sample}/${sample}.dedup_GC.sorted.bam ${data_dir}/${sample}/${sample}.dedup_GC.bam
    rm ${data_dir}/${sample}/${sample}.dedup_GC.bam
    samtools index ${data_dir}/${sample}/${sample}.dedup_GC.sorted.bam

done