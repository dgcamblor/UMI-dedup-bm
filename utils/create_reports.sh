#!/bin/bash

#-------------------------------------------------------------------------------
# create_reports.sh
#-------------------------------------------------------------------------------
# Create IGV reports to evaluate the outliers in the UMIvar evaluation
#-------------------------------------------------------------------------------

create_report results/umivar_test/umivar_outliers.bed \
    --genome hg38 \
    --flanking 100 \
    --tracks \
        results/umivar_test/umivar_outliers.bed \
        data/test_umivar/EP7/EP7_UV.sorted.bam \
        data/test_umivar/EP9/EP9_UV.sorted.bam \
        data/test_umivar/EP11/EP11_UV.sorted.bam \
        data/test_umivar/EP13/EP13_UV.sorted.bam \
        data/test_umivar/EP49/EP49_UV.sorted.bam \
    --output results/umivar_test/umivar_outliers_report.html
