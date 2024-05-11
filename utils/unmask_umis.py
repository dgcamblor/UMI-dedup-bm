#-------------------------------------------------------------------------------
# unmask_umis.py
#-------------------------------------------------------------------------------
# This script takes as input two FASTQ files (paired end reads) and outputs the
# same files with unmasked UMIs. A dictionary is used to store the UMI masked 
# sequence and a 12 bp UMI is randomly generated for that sequence.
#-------------------------------------------------------------------------------

import sys
import random
import argparse
import pickle
import subprocess

masked_to_umi = {}


def random_umi():
    umi = ""
    for i in range(12):
        umi += random.choice("ATCG")
    return umi


def process_fastq(fastq1, fastq2, output1, output2, n_reads=None):
    # Open the files
    f1 = open(fastq1, "r")
    f2 = open(fastq2, "r")
    o1 = open(output1, "w")
    o2 = open(output2, "w")

    processed_reads = 0

    # Read the files line by line
    for line1, line2 in zip(f1, f2):
        # Process the first line
        if line1.startswith("@"):
            # Split the line by " "
            fields1 = line1.split(" ")
            fields2 = line2.split(" ")

            if fields1[0] != fields2[0]:
                print("ERROR: The FASTQ files are not paired end reads.")
                sys.exit(42)

            # Remove the third field
            fields1.pop(2)
            fields2.pop(2)
            # Join the fields with ":"
            fields1 = ":".join(fields1)
            fields2 = ":".join(fields2)

            masked_umi = fields1.split(":")[-1]
            if masked_umi in masked_to_umi:
                umi = masked_to_umi[masked_umi]
            else:
                umi = random_umi()
                masked_to_umi[masked_umi] = umi

            fields1 = f"{fields1}_{umi}\n"
            fields2 = f"{fields2}_{umi}\n"

            # Write the processed lines
            o1.write(fields1)
            o2.write(fields2)

            if n_reads:
                processed_reads += 1
                sys.stdout.write(f"\rProcessed {processed_reads} reads ({processed_reads / n_reads * 100:.2f}%)")
                sys.stdout.flush()

        else:
            o1.write(line1)
            o2.write(line2)

    # Close the files
    f1.close()
    f2.close()
    o1.close()
    o2.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastq1", help="FASTQ file 1", required=True)
    parser.add_argument("--fastq2", help="FASTQ file 2", required=True)
    parser.add_argument("--output1", help="Output FASTQ file 1", required=True)
    parser.add_argument("--output2", help="Output FASTQ file 2", required=True)
    args = parser.parse_args()

    sra_id = args.fastq1.split("_")[0].split("/")[-1]

    print(f"Initiated process_fastq.py for {sra_id}...")

    # Get the number of reads in the FASTQ file
    cmd = "awk 'END {print NR/4}' " + args.fastq1

    n_reads = int(subprocess.check_output(cmd, shell=True))

    process_fastq(args.fastq1, args.fastq2, args.output1, args.output2, n_reads)

    # Save the dictionary of masked to unmasked UMIs
    with open(f"{sra_id}.pkl", "wb") as f:
        pickle.dump(masked_to_umi, f)


if __name__ == "__main__":
    main()