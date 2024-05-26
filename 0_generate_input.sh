source config.conf
eval "$(conda shell.bash hook)"; conda activate dedup_bm
shopt -s expand_aliases

#-------------------------------------------------------------------------------
# Process BRP true variants
#-------------------------------------------------------------------------------

# Create on-target true variants (SEQC2 known positives and BRP sequenced regions)
echo -e "Creating on-target true variants for BRP"
bedtools intersect -a metadata/KnownPositives_hg19.vcf.gz -b metadata/LBx_BRP_hg19.bed -header > metadata/KnownPositives_hg19_target.vcf

# Normalize VCF
echo -e "Normalizing the intersection"
bcftools norm -m-any --check-ref w -f $REF_GENOME_HG19 \
    -o metadata/KnownPositives_hg19_target_norm.vcf \
    metadata/KnownPositives_hg19_target.vcf

#-------------------------------------------------------------------------------
# Create input for UMIvar
#-------------------------------------------------------------------------------

# UMIvar target: exonic regions of cancer-related genes

echo -e "Creating UMIvar target regions"
./utils/create_umivar_target.sh \
    "/nfs/home/panel_designs/HyperExome/target.bed" \
    "/media/scratchDELLEMC/dedup_bm/metadata/cancerGeneList.tsv" \
    "metadata/umivar_target.bed" \
    100

echo -e "Intersecting dbSNP v155 (oncoKB/COSMIC annotations) to target regions for UMIvar"
bedtools intersect \
    -a metadata/dbSNP155.cosmic.oncokb.annotation.vcf.gz \
    -b metadata/umivar_target.bed \
    -header > metadata/dbsnp_v155_onco_annot.ontarget.vcf

echo -e "Normalizing the intersection"
bcftools_nfs norm -m-any --check-ref w -f $REF_GENOME_GRCH38 -o metadata/dbsnp_v155_onco_annot.ontarget.norm.vcf metadata/dbsnp_v155_onco_annot.ontarget.vcf

# UMIvar testing variants: random variants from dbSNP v155

echo -e "Selecting variants for UMIvar testing"
python3 utils/select_vars_umivar.py \
    --vcf_in "metadata/dbsnp_v155_onco_annot.ontarget.norm.vcf" \
    --csv_out "metadata/umivar_test_input.csv" \
    --n_vars 1320 \
    --min_af 0.01 \
    --max_af 1.00

echo -e "Filtering for suitable variants"
python3 utils/filter_variants.py \
    --csv_in "metadata/umivar_test_input.csv" \
    --csv_out "metadata/umivar_test_input_filt.csv"

echo -e "Final number of variants: $(wc -l metadata/umivar_test_input_filt.csv)"

mv metadata/umivar_test_input.csv metadata/umivar_test_input_unfiltered.csv
mv metadata/umivar_test_input_filt.csv metadata/umivar_test_input.csv

# UMIvar benchmarking dataset: oncogenic / likely oncogenic variants

echo -e "Selecting oncogenic / likely oncogenic variants"
bcftools view -i 'INFO/Oncogenic="Likely_oncogenic" || INFO/Oncogenic="Oncogenic"' metadata/dbsnp_v155_onco_annot.ontarget.norm.vcf > metadata/dbsnp_v155_onco_annot.ontarget.norm.oncogenic.vcf

echo -e "Selecting variants for UMIvar dataset"
python3 utils/select_vars_umivar.py \
    --vcf_in "metadata/dbsnp_v155_onco_annot.ontarget.norm.oncogenic.vcf" \
    --csv_out "metadata/umivar_input.csv" \
    --n_vars 1320 \
    --min_af 0.01 \
    --max_af 0.10

echo -e "Filtering for suitable variants"
python3 utils/filter_variants.py \
    --csv_in "metadata/umivar_input.csv" \
    --csv_out "metadata/umivar_input_filt.csv"

echo -e "Final number of variants: $(wc -l metadata/umivar_input_filt.csv)"

mv metadata/umivar_input.csv metadata/umivar_input_unfiltered.csv
mv metadata/umivar_input_filt.csv metadata/umivar_input.csv