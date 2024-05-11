#' Extract variants from a VCF file
#' This function is used to extract the 230 true variants for BRP, and also
#' the variants called with each different pipeline.
#' 
#' @param variant_file Path to the VCF file
#' @param af_meta The meta data where the allele frequency is stored. Can be either "gt" or "fix"
#' @param af_name The name of the allele frequency field in the meta data ("AF", "VAF"...)
#' @return A list of GeneticVariant objects

source("R/GeneticVariant.R")

extract_variants <- function(variant_file, af_meta="gt", af_name="AF") {
  vcf <- read.vcfR(variant_file, verbose = FALSE)
  vcf_fix <- vcf@fix %>% as.data.frame()
  vcf_gt <- vcf@gt %>% as.data.frame()
    
  chroms <- vcf_fix %>% pull(CHROM)
  pos <- vcf_fix %>% pull(POS)
  refs <- vcf_fix %>% pull(REF)
  alts <- vcf_fix %>% pull(ALT)
  pass <- vcf_fix %>% pull(FILTER)
  
  # Get the allele frequency of the mutation
  if (af_meta == "gt") {
    afs <- sapply(str_split(vcf_gt[,2], ":"), function(x){x[[which(str_detect(str_split(vcf_gt[1,1], ":")[[1]], af_name))]]})
  } 
  
  else if (af_meta == "fix") {
    vcf_info <- vcf_fix %>% pull(INFO)
    afs <- sapply(str_split(vcf_info, ";"), function(x){x[[which(str_detect(str_split(x, ";")[[1]], af_name))]]})
  }
  
  # Remove everything that is not a number or a period
  afs <- gsub("[^0-9.]", "", afs)
  
  # Create GeneticVariant objects
  variants <- sapply(1:length(chroms), function(i) {
    new("GeneticVariant", CHR = chroms[i], POS = as.numeric(pos[i]), REF = refs[i], ALT = alts[i], VAF = as.numeric(afs[i]))
  })
  
  # Filter out variants that do not pass (FILTER column)
  variants <- variants[pass == "PASS" | is.na(pass)]
  
  print(paste0("Extracted ", length(variants), " variants from ", variant_file))
  print(paste0("Filtered out ", length(chroms) - length(variants), " variants that did not PASS"))
  
  return(variants)
}