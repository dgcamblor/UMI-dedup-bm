# Pipeline to benchmark UMI deduplication tools in variant calling

## Introduction

This pipeline benchmarks the effect in the downstream variant calling results of performing UMI-based read deduplication with different software: Picard MarkDuplicates (using `--BARCODE_TAG`), UMI-tools, fgbio, UMICollapse and gencore. The benchmarking is performed in two datasets:

- A real dataset from deep sequencing (BRP data) belonging to the [SEQC2 Liquid Biopsy study](https://www.nature.com/articles/s41597-022-01276-8).
- A simulated dataset created with UMIvar.

## Description

A brief explanation of the scripts used in the pipeline is provided below:

- `launch_setup.sh`: Downloads the necessary tools for the pipeline execution.
- `launch_pipeline.sh`: Launches the pipeline (numbered scripts below).

- `0_generate_input.sh`: Generates the random variants to be used in the simulated dataset, and processes the variants of the real dataset.
- `1_test_umivar.sh`: Simulates a testing dataset with UMIvar to study the performance of this tool, and performs variant calling with VarDict. The `launch_comp.sh` script is used to create a comparative table of the results.
- `2_launch_umivar.sh`: Creates a simulated dataset with UMIvar for further use in the benchmarking.
- `3_download_and_process.sh`: Downloads the real dataset from the BRP sequencing and processes it. The UMIs in the data are masked, so an additional step is required to unmask them. 
- `4_dedup.sh`: Performs the UMI deduplication with the tools mentioned above.
- `5_variant_calling.sh`: Calls variants with Mutect2, Lofreq and VarDict. Minimal preprocessing is performed before the variant calling (`BQSR`).
- `6_gather_stats.sh`: Gathers additional statistics of the samples.
- `7_benchmark.qmd`: A Quarto Markdown file to analyze the results of the benchmarking.

Additionally, the `utils` folder contains several miscellaneous scripts used in the pipeline, and the `R` folder contains the functions and classes used in the `7_benchmark.qmd` script.