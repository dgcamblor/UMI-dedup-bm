#-------------------------------------------------------------------------------
# umi2seq.py
#-------------------------------------------------------------------------------
# This scripts adds the UMI information stored in a header to the sequence.
# Qualities are set to "I".
#-------------------------------------------------------------------------------

import argparse

QUALITY = "I"

parser = argparse.ArgumentParser(description='Extracts UMI in headers from a FASTQ file and adds it to the sequence.')
parser.add_argument("-f", "--file", help="Input file in FASTQ format", required=True)
args = parser.parse_args()


out_path = "".join(args.file.split(".")[0]) + "_UIS.fastq"

with(
    open(args.file, 'r') as f_in,
    open(out_path, 'w') as f_out
):
    for i, line in enumerate(f_in):
        if i % 4 == 0:
            header = line.strip()
            umi = header.split("_")[-1]

        elif i % 4 == 1:
            seq = line.strip()

        elif i % 4 == 2:
            plus = line.strip()

        elif i % 4 == 3:
            qual = line.strip()

            f_out.write(header + "\n")
            f_out.write(umi + seq + "\n")
            f_out.write(plus + "\n")
            f_out.write(QUALITY * len(umi) + qual + "\n")