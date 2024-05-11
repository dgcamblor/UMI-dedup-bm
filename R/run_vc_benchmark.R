#' Run a variant calling benchmark in the deduplication benchmarking pipeline
#' 
#' @param variant_caller The variant caller to be benchmarked (to extract the results)
#' @param profile The profile to be used for the benchmarking (if BRP, known false positives are used)
#' 
#' @return A list containing the benchmark results and statistics

run_vc_benchmark <- function(variant_caller, profile) {
  benchmarks <- list()
  stats <- list()

  for (accession in accessions) {
    # List the files containing the results from the variant calling
    file_list <- list.files(paste0(results_dir, "/", accession, "/", variant_caller), pattern = ".*target.*_norm.vcf*", full.names = TRUE)
    file_list <- file_list[!grepl(".tbi", file_list)]  # Remove index files

    # For each file listed, extract the variants to a GeneticVariant list    
    called_variants <- lapply(file_list, extract_variants); names(called_variants) <- basename(file_list)
    
    # Logic to determine false positives: known or all non-true positives
    if (profile == "BRP") {
      # Use bedtools intersect to extracts variants called in known false negative regions
      for (file in file_list) {
        system(paste("bedtools intersect -a", file, "-b metadata/KnownNegatives_hg19.bed.gz -header >", gsub(".vcf.gz", "_known_false_positives.vcf", file)))
      }
      
      # Extract the variants from the known false positive regions to a GeneticVariant list
      known_false_positives <- lapply(gsub(".vcf.gz", "_known_false_positives.vcf", file_list), extract_variants)
      names(known_false_positives) <- basename(file_list)

      benchmark_results <- vc_benchmark(true_variants, called_variants, known_false_positives)
    } else if (profile == "UMIVAR") {
      # Use all non-true positives as false positives
      benchmark_results <- vc_benchmark(true_variants, called_variants)
    } else {
      stop("Invalid profile. Please provide a valid profile.")
    }
    
    # Store the benchmark results in the benchmarks list
    benchmarks[[accession]] <- benchmark_results
    
    # Store the statistics in the stats list
    df <- as.data.frame(benchmark_results)
    stats[[accession]] <- df
  }

  rm(benchmark_results)
  rm(df)

  # Combine all dataframes into a single dataframe
  stats <- do.call(rbind, stats)
  return(list(benchmarks = benchmarks, stats = stats))
}