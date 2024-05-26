#!/bin/bash

#-------------------------------------------------------------------------------
# launch_pipeline.sh
#-------------------------------------------------------------------------------
# Execute the steps of the benchmarking pipeline.
#-------------------------------------------------------------------------------

# launch_setup.sh
./0_generate_input.sh

./1_test_umivar.sh
./launch_comp.sh

./2_launch_umivar.sh
./3_download_and_process.sh
./4_dedup.sh
./5_variant_calling.sh
./6_gather_stats.sh