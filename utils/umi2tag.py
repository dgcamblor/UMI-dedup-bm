#-------------------------------------------------------------------------------
# umi2tag.py
#-------------------------------------------------------------------------------
# This script takes a BAM as an input, in which read names have the UMI appended 
# after "_", and adds a RX tag to the BAM file with the UMI sequence. 
# The output is a BAM file with the RX tag added.
#-------------------------------------------------------------------------------

import pysam
import argparse

def umi2tag(in_bam, out_bam):
    in_bam = pysam.AlignmentFile(in_bam, "rb")
    out_bam = pysam.AlignmentFile(out_bam, "wb", template=in_bam)

    for read in in_bam.fetch(until_eof=True):
        umi = read.query_name.split("_")[-1]
        read.set_tag("RX", umi)
        out_bam.write(read)

    in_bam.close()
    out_bam.close()

def main():
    print(f"Adding UMI to tags in {args.input}...") 

    umi2tag(args.input, args.output)
    
    print("Done!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extracts UMI in headers from a BAM file and adds it to the sequence.')
    parser.add_argument("-i", "--input", help="Input file in BAM format", required=True)
    parser.add_argument("-o", "--output", help="Output file in BAM format", required=True)
    args = parser.parse_args()

    main()