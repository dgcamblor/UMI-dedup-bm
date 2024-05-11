#!/usr/bin/env python3

#-------------------------------------------------------------------------------
# create_umivar_input.py
#-------------------------------------------------------------------------------
# From a VCF file, create a CSV file in a format suitable for UMIvar input.
# A random allele frequency between 0.01 and 0.10 is assigned to each variant.
#-------------------------------------------------------------------------------

import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Create a CSV file with random allele frequencies for UMIvar input.")
parser.add_argument("--vcf_in", help="Path to the VCF file.")
parser.add_argument("--csv_out", help="Path to the CSV file (UMIvar input) to be created.")
parser.add_argument("--min_af", help="Minimum allele frequency.", type=float)
parser.add_argument("--max_af", help="Maximum allele frequency.", type=float)
args = parser.parse_args()


np.random.seed(12)

allele_freqs = np.random.uniform(args.min_af, args.max_af, 100)

# Read the VCF file
with open(args.vcf_in, "r") as in_file:
    with open (args.csv_out, "w") as out_file:
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

                out_file.write(f"{chrom},{pos},{ref},{alt},{af}\n")