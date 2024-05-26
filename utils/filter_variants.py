#!/usr/bin/env python3

#-------------------------------------------------------------------------------
# filter_variants.py
#-------------------------------------------------------------------------------
# This script removes variants that are located within 10 bp of each other.
# Overlaps are avoided for adequate UMIvar evaluation.
#-------------------------------------------------------------------------------

import argparse

parser = argparse.ArgumentParser(description="Filter variants that are located within 10 bp of each other.")
parser.add_argument("--csv_in", help="Path to the input CSV file.")
parser.add_argument("--csv_out", help="Path to the output CSV file.")
args = parser.parse_args()

with open(args.csv_in, 'r') as f:
    lines = f.readlines()

with open(args.csv_out, 'w') as f:
    chrom_prev = ""
    pos_prev = 0
    for line in lines:
        chrom, pos, ref, alt, score = line.strip().split(',')
        pos = int(pos)
        if chrom != chrom_prev:
            f.write(f'{chrom},{pos},{ref},{alt},{score}\n')
            chrom_prev = chrom
            pos_prev = pos
        elif chrom == chrom_prev and abs(pos - pos_prev) > 10:
            f.write(f'{chrom},{pos},{ref},{alt},{score}\n')
            chrom_prev = chrom
            pos_prev = pos
        else:
            chrom_prev = chrom
            pos_prev = pos