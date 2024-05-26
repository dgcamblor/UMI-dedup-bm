#!/usr/bin/env python3

#-------------------------------------------------------------------------------
# select_vars_umivar.py
#-------------------------------------------------------------------------------
# From a VCF file, select random variants to create a CSV file in a format 
# suitable for UMIvar input. Allele frequencies are randomly generated.
#-------------------------------------------------------------------------------

import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Create a CSV file with random allele frequencies for UMIvar input.")
parser.add_argument("--vcf_in", help="Path to the VCF file.")
parser.add_argument("--csv_out", help="Path to the CSV file (UMIvar input) to be created.")
parser.add_argument("--n_vars", help="Number of variants to select.", type=int)
parser.add_argument("--min_af", help="Minimum allele frequency.", type=float)
parser.add_argument("--max_af", help="Maximum allele frequency.", type=float)
args = parser.parse_args()


np.random.seed(12)

allele_freqs = np.random.uniform(args.min_af, args.max_af, 100)

variants = []

# Read the VCF file
with open(args.vcf_in, "r") as in_file:
    for line in in_file:
        if line.startswith("#"):
            continue
        else:
            line = line.strip().split("\t")
            chrom = line[0]
            pos = line[1]
            ref = line[3]
            alt = line[4]
            af = np.random.choice(allele_freqs)
            variants.append((chrom, pos, ref, alt, af))

# Filter variants suitable for evaluation
variants = [var for var in variants if len(var[2]) <= 10 and len(var[3]) <= 10]
variants = [var for var in variants if len(var[2]) == 1 or len(var[3]) == 1]
variants = [var for var in variants if var[2][0] == var[3][0] or (len(var[2]) == 1 and len(var[3]) == 1)]

# Select random variants
selected_i = np.random.choice(len(variants), args.n_vars, replace=False)
selected_vars = [variants[i] for i in selected_i]

# Separate chrX variants
chrX_vars = [var for var in selected_vars if var[0] == "chrX"]
selected_vars = [var for var in selected_vars if var[0] != "chrX"]

# Order variants by chromosome and position
selected_vars.sort(key=lambda x: (int(x[0][3:]), int(x[1])))
chrX_vars.sort(key=lambda x: int(x[1]))

# Add back the chrX variants
selected_vars += chrX_vars

# Write the selected variants to the output file
with open(args.csv_out, "w") as out_file:
    for var in selected_vars:
        out_file.write(f"{var[0]},{var[1]},{var[2]},{var[3]},{var[4]}\n")