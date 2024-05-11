#' Extract variants from the UMIvar input
#' 
#' @param variant_file File containing the variants
#' @return List of GeneticVariant objects

source("R/GeneticVariant.R")

extract_variants_umivar <- function(variant_file){
  var_file <- read.csv(variant_file)
  
  variants <- list()
  
  for (i in 1:nrow(var_file)){
    # Extract the variant information
    variants <- append(variants, new("GeneticVariant", CHR = var_file[i, 1], POS = var_file[i, 2], REF = var_file[i, 3], ALT = var_file[i, 4], VAF = var_file[i, 5]))
  }
  
  return(variants)
}