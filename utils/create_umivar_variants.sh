#!/bin/bash

shopt -s expand_aliases

create_umivar_variants() {
    local genome_file="$1"  #$ The reference genome file
    local input_vcf="$2"  #$ The input VCF file (a VCF file from a variant database)
    local output_vcf="$3"  #$ The output VCF file (a VCF file with the variants to be used for UMIvar)
    local target_bed="$4"  #$ A BED file containing the target regions)
    local random_fraction="$5"  #$ The fraction of variants to be selected randomly

    dir_vcf=$(dirname "$input_vcf")
    name_vcf=$(basename "$input_vcf")
    name_vcf=$(echo "$name_vcf" | cut -d'.' -f1)

    # If the input VCF does not contain .norm, normalize it with bcftools norm
    if [[ "$input_vcf" != *".norm."* ]]; then
        echo "Normalizing..."
        bcftools norm -m-any -f "$genome_file" -o "$dir_vcf/$name_vcf.norm.vcf.gz" "$input_vcf"
        input_vcf="$dir_vcf/$name_vcf.norm.vcf.gz"
    
    else 
        echo "The input VCF file is already normalized."

    fi

    # If the input VCF does not contain .ontarget, intersect it with the UMIvar target
    if [[ "$input_vcf" != *".ontarget."* ]]; then
        echo "Intersecting with UMIvar target..."
        bedtools intersect -a "$input_vcf" -b "$target_bed" -header > "$dir_vcf/$name_vcf.ontarget.norm.vcf"
        input_vcf="$dir_vcf/$name_vcf.ontarget.norm.vcf"

    else 
        echo "The input VCF file is already intersected with the UMIvar target."

    fi

    # Downsampling the variants
    echo "Downsampling the variants..."
    gatk SelectVariants -V "$input_vcf" -O "$output_vcf" --select-random-fraction "$random_fraction"
}

create_umivar_variants "$1" "$2" "$3" "$4" "$5"