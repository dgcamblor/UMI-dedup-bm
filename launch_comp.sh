#!/bin/bash

#-------------------------------------------------------------------------------
# launch_comp.sh
#-------------------------------------------------------------------------------
# Compare the results obtained from the UMIvar testing (VarDict calling in test
# variants) vs. the input of the user.
#-------------------------------------------------------------------------------

source config.conf

# Create the conda environment
eval "$(conda shell.bash hook)"; conda activate dedup_bm

SAMPLES_UMIVAR_TEST_STR=$(IFS=,; echo "${SAMPLES_UMIVAR_TEST[*]}")  # Join array elements with comma

python3 utils/comp_umivar_test.py \
    --umivar_input "/media/scratchDELLEMC/dedup_bm/metadata/umivar_test_input.csv" \
    --umivar_dir "/media/scratchDELLEMC/dedup_bm/data/umivar_test" \
    --vardict_dir "/media/scratchDELLEMC/dedup_bm/results/umivar_test" \
    --samples $SAMPLES_UMIVAR_TEST_STR \
    --output "/media/scratchDELLEMC/dedup_bm/results/umivar_test/test_umivar_results.csv"