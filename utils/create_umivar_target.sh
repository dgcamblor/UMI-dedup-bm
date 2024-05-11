#!/bin/bash

#-------------------------------------------------------------------------------
# create_umivar_target.sh
#-------------------------------------------------------------------------------

function create_umivar_target() {
    local exome_bed_file="$1"  #$ A bed file containing the exome regions
    local cancer_genes_file="$2"  #$ A tsv file containing cancer genes
    local extracted_regions_file="$3"  #$ Output file to store the extracted regions
    local min_len="$4"  #$ Minimum length of the extracted regions

    # Get cancer genes
    awk '$10=="Yes" && $16=="Yes" {print $1}' "$cancer_genes_file" | sort > "cancer_gene_list.txt"

    # Filter target bed file based on unique genes
    grep -wf "cancer_gene_list.txt" "$exome_bed_file" > "$extracted_regions_file"

    rm "cancer_gene_list.txt"

    # Filter the extracted regions based on minimum length
    local extrated_regions_filtered_file="${extracted_regions_file%.bed}_filtered.bed"
    awk -v min_len="$min_len" '{if ($3-$2 > min_len) print $0}' "$extracted_regions_file" > "$extrated_regions_filtered_file"

    mv "$extrated_regions_filtered_file" "$extracted_regions_file"

    echo "Extraction complete. Results are in $extracted_regions_file"
}

# Call the function with arguments
create_umivar_target "$1" "$2" "$3" "$4"
