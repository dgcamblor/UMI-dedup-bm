#!/usr/bin/env python3

#-------------------------------------------------------------------------------
# comp_umivar_test.py
#-------------------------------------------------------------------------------
# Compare the VarDict variant calling results with the expected results (input)
# for the UMIvar test variants.
#-------------------------------------------------------------------------------

import argparse

parser = argparse.ArgumentParser(description="Compare the VarDict variant calling results with the expected results for the UMIvar test variants.")
parser.add_argument("--umivar_input", help="Path to the CSV file with the UMIvar variants.")
parser.add_argument("--samples", help="List of samples to be tested.")
parser.add_argument("--umivar_dir", help="Base path to the directory with the UMIvar CSV files.")
parser.add_argument("--vardict_dir", help="Base path to the directory with the VarDict VCF files.")
parser.add_argument("--output", help="Path to the output file.")
args = parser.parse_args()

# Extract the list of variants to test
test_vars_results = {}

with open(args.umivar_input, "r") as f:
    for line in f:
        # Separate the line into a list of strings
        line = line.strip().split(",")

        chr = line[0]
        pos = line[1]
        ref = line[2]
        alt = line[3]


        id = chr + ":" + pos + ref + ">" + alt  # Define a variant identifier
        af = line[4]

        test_vars_results[id] = [af]

samples = args.samples.split(",")

# Add the UMIvar estimates for each sample
for sample in samples:
    # Extract the list of variants called by UMIvar
    umivar_results = {}
    with open(args.umivar_dir + "/" + sample + "/" + sample + "_UV.csv", "r") as f:
        # Skip the header
        f.readline()

        for line in f:
            # Separate the line into a list of strings
            line = line.strip().split(",")

            chr = line[0]
            pos = line[1]
            ref = line[2]
            alt = line[3]

            id = chr + ":" + pos + ref + ">" + alt  # Again retrieve the variant identifier
            af = line[4]

            umivar_results[id] = af

    # Check for each variant whether it has been called by UMIvar
    for var in test_vars_results:
        if var in umivar_results:
            # Append the AF to the list
            test_vars_results[var].append(umivar_results[var])
        else:
            test_vars_results[var].append("0")  # Codify non-called variants with AF 0

# Check for each sample whether the test variants have been called by VarDict
for sample in samples:
    # Extract the list of variants called by VarDict
    vardict_results = {}
    with open(args.vardict_dir + "/" + sample + "_vardict_norm.vcf", "r") as f:
        for line in f:
            # Skip the header
            if line[0] == "#":
                continue

            # Separate the line into a list of strings
            line = line.strip().split("\t")

            chr = line[0]
            pos = line[1]
            ref = line[3]
            alt = line[4]

            id = chr + ":" + pos + ref + ">" + alt  # Again retrieve the variant identifier
            af = line[7].split(";")[4].split("=")[1]

            vardict_results[id] = af

    # Check for each variant whether it has been called by VarDict
    for var in test_vars_results:
        if var in vardict_results:
            # Append the AF to the list
            test_vars_results[var].append(vardict_results[var])
        else:
            test_vars_results[var].append("0")  # Codify non-called variants with AF 0

samples = [sample + "_umivar" for sample in samples]
samples += [sample + "_vardict" for sample in args.samples.split(",")]

# Write the results to a file
with open(args.output, "w") as f:
    f.write("var,af," + ",".join(samples) + "\n")
    for var in test_vars_results:
        f.write(var + "," + ",".join(test_vars_results[var]) + "\n")