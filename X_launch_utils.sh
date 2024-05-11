source config.conf
eval "$(conda shell.bash hook)"; conda activate dedup_bm

./utils/create_umivar_target.sh \
    "/nfs/home/panel_designs/HyperExome/target.bed" \
    "/media/scratchDELLEMC/dedup_bm/metadata/cancerGeneList.tsv" \
    "metadata/umivar_target.bed" \
    100

#-------------------------------------------------------------------------------
# UMIvar benchmarking
#-------------------------------------------------------------------------------
./utils/create_umivar_variants.sh \
    "/nfs/home/references/genomes/human/GRCh38/Homo_sapiens.GRCh38.fa" \
    "metadata/dbsnp_v156_20231113.ontarget.norm.vcf" \
    "metadata/umivar_input.vcf" \
    "metadata/umivar_target.bed" \
    0.0006

python3 utils/create_umivar_input.py \
    --vcf_in "/media/scratchDELLEMC/dedup_bm/metadata/umivar_input.vcf" \
    --csv_out "/media/scratchDELLEMC/dedup_bm/metadata/umivar_input.csv" \
    --min_af 0.01 \
    --max_af 0.10  # 0.10

#-------------------------------------------------------------------------------
# UMIvar test
#-------------------------------------------------------------------------------
./utils/create_umivar_variants.sh \
    "/nfs/home/references/genomes/human/GRCh38/Homo_sapiens.GRCh38.fa" \
    "metadata/dbsnp_v156_20231113.ontarget.norm.vcf" \
    "metadata/umivar_test_input.vcf" \
    "metadata/umivar_target.bed" \
    0.001

python3 utils/create_umivar_input.py \
    --vcf_in "/media/scratchDELLEMC/dedup_bm/metadata/umivar_test_input.vcf" \
    --csv_out "/media/scratchDELLEMC/dedup_bm/metadata/umivar_test_input.csv" \
    --min_af 0.01 \
    --max_af 1.0  # 0.10

SAMPLES_UMIVAR_TEST_STR=$(IFS=,; echo "${SAMPLES_UMIVAR_TEST[*]}")

python3 utils/comp_umivar_test.py \
    --umivar_input "/media/scratchDELLEMC/dedup_bm/metadata/umivar_test_input.csv" \
    --samples $SAMPLES_UMIVAR_TEST_STR \
    --vardict_dir "/media/scratchDELLEMC/dedup_bm/data/test_umivar/" \
    --output "/media/scratchDELLEMC/dedup_bm/test_umivar_results.csv"
