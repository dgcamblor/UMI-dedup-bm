#' Execute a variant calling benchmark with a set of true variants and called variants
#' 
#' @param true_variants A list of true variants to be used as a reference (list of GeneticVariant objects)
#' @param called_variants A list of pipelines with their corresponding list of called variants (list of lists of GeneticVariant objects)
#' 
#' @return A BenchmarkResults object that contains the benchmarking statistics

source("R/BenchmarkResults.R")

vc_benchmark <- function(true_variants, called_variants, known_false_positives=NULL) {
  sample_name <- str_split(names(called_variants)[1], "_")[[1]][1]
  pipeline_names <- sapply(str_split(names(called_variants), "_"), function(x) x[2])
  
  global_tp <- list()
  global_fn <- list()
  global_fp <- list()
  
  for (file in names(called_variants)) {
    pipeline_name <- str_split(file, "_")[[1]][2]
    
    file_tp <- list()
    file_fn <- list()
    file_fp <- list()
    
    # Compute the true positives and false negatives
    for (true_variant in true_variants) {
      found <- FALSE  # Set to TRUE if the variant is found in the called variants
      
      for (called_variant in called_variants[[file]]) {
        if (true_variant == called_variant) {
          file_tp <- append(file_tp, true_variant)
          found <- TRUE
        }
      }
      
      # If the variant is not found, it is a false negative
      if (!found) {
        file_fn <- append(file_fn, true_variant)
      }
    }
      
    # Compute the false positives: if no negative regions are provided, we assume that all non-tp variants are fp
    if (is.null(known_false_positives)){
      for (called_variant in called_variants[[file]]) {
        found <- FALSE  # Set to TRUE if the variant is found in the true variants
        
        for (true_variant in true_variants) {
          if (true_variant == called_variant) {
            found <- TRUE
          }
        }
        
        # If the variant is not found, it is a false positive
        if (!found) {
          file_fp <- append(file_fp, called_variant)
        }
      }
    }
    else{
      # If negative regions are provided, sort the negative regions by position
      file_fp <- known_false_positives[[file]]
    }
    
    # Append the results to the global variables
    global_tp[[pipeline_name]] <- file_tp
    global_fn[[pipeline_name]] <- file_fn
    global_fp[[pipeline_name]] <- file_fp
      
  }
  
  # Create the object
  benchmark_results <- new("BenchmarkResults", 
                           sample_name = sample_name, 
                           pipeline_names = pipeline_names, 
                           true_positives = global_tp, 
                           false_positives = global_fp, 
                           false_negatives = global_fn)
  
  return(benchmark_results)
}