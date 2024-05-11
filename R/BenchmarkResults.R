#' Class for working with benchmark results from variant calling using different pipelines
#' 
#' @slot sample_name The name of the sample
#' @slot pipeline_names The names of the pipelines that are being compared
#' @slot true_positives A list of true positives for each pipeline
#' @slot false_positives A list of false positives for each pipeline
#' @slot false_negatives A list of false negatives for each pipeline

setClass(
  "BenchmarkResults",
  representation(
    sample_name = "character",
    pipeline_names = "character",
    true_positives = "list",
    false_positives = "list",
    false_negatives = "list"
  )
)

# Convert to data frame (with counts)
setMethod("as.data.frame", "BenchmarkResults", function(x) {
  tp_count <- get_tp_count(x)
  fp_count <- get_fp_count(x)
  fn_count <- get_fn_count(x)
  
  df <- data.frame(
    sample_name = x@sample_name,
    pipeline_names = x@pipeline_names,
    tp_count = tp_count,
    fp_count = fp_count,
    fn_count = fn_count,
    recall = get_recall(x),
    precision = get_precision(x),
    f1 = get_f1(x),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  return(df)
})

setGeneric("get_exclusive_variants", function(object, variant_type, target_pipeline, comparison_pipelines) {
  standardGeneric("get_exclusive_variants")
})

setMethod("get_exclusive_variants", "BenchmarkResults", function(object, variant_type, target_pipeline, comparison_pipelines) {
  if (!(variant_type %in% c("tp", "fp", "fn"))) {
    stop("Invalid variant_type. It should be one of 'tp', 'fp', or 'fn'")
  }
  
  if (!(target_pipeline %in% object@pipeline_names)) {
    stop("Target pipeline not found in the benchmark results.")
  }
  
  if (length(comparison_pipelines) == 0 || any(!(comparison_pipelines %in% object@pipeline_names))) {
    stop("Invalid comparison pipelines. Please provide valid pipeline names.")
  }
  
  # Get the index of the target and comparison pipelines
  target_index <- which(object@pipeline_names == target_pipeline)
  comparison_indices <- which(object@pipeline_names %in% comparison_pipelines)
  
  # Retrieve exclusive variants
  exclusive_variants <- list()
  for (i in comparison_indices) {
    if (variant_type == "tp") {
      exclusive_variants[[object@pipeline_names[i]]] <- setdiff(object@true_positives[[target_index]], object@true_positives[[i]])
    } else if (variant_type == "fp") {
      exclusive_variants[[object@pipeline_names[i]]] <- setdiff(object@false_positives[[target_index]], object@false_positives[[i]])
    } else if (variant_type == "fn") {
      exclusive_variants[[object@pipeline_names[i]]] <- setdiff(object@false_negatives[[target_index]], object@false_negatives[[i]])
    }
  }
  
  return(exclusive_variants)
})

# Getters ----------------------------------------------------------------------
setGeneric("get_sample_name", function(object) standardGeneric("get_sample_name"))
setGeneric("get_pipeline_names", function(object) standardGeneric("get_pipeline_names"))

setGeneric("get_tp", function(object) standardGeneric("get_tp"))
setGeneric("get_tp_count", function(object) standardGeneric("get_tp_count"))
setGeneric("get_fp", function(object) standardGeneric("get_fp"))
setGeneric("get_fp_count", function(object) standardGeneric("get_fp_count"))
setGeneric("get_fn", function(object) standardGeneric("get_fn"))
setGeneric("get_fn_count", function(object) standardGeneric("get_fn_count"))

setGeneric("get_recall", function(object) standardGeneric("get_recall"))
setGeneric("get_precision", function(object) standardGeneric("get_precision"))
setGeneric("get_f1", function(object) standardGeneric("get_f1"))

setMethod("get_sample_name", "BenchmarkResults", function(object) {
  return(object@sample_name)
})

setMethod("get_pipeline_names", "BenchmarkResults", function(object) {
  return(object@pipeline_names)
})

setMethod("get_tp", "BenchmarkResults", function(object) {
  return(object@true_positives)
})

setMethod("get_tp_count", "BenchmarkResults", function(object) {
  return(sapply(object@true_positives, length))
})

setMethod("get_fp", "BenchmarkResults", function(object) {
  return(object@false_positives)
})

setMethod("get_fp_count", "BenchmarkResults", function(object) {
  return(sapply(object@false_positives, length))
})

setMethod("get_fn", "BenchmarkResults", function(object) {
  return(object@false_negatives)
})

setMethod("get_fn_count", "BenchmarkResults", function(object) {
  return(sapply(object@false_negatives, length))
})

setMethod("get_recall", "BenchmarkResults", function(object) {
  recalls <- sapply(1:length(object@pipeline_names), function(i) {
    return(length(object@true_positives[[i]]) / (length(object@true_positives[[i]]) + length(object@false_negatives[[i]])))
    })
  names(recalls) <- object@pipeline_names
  return(recalls)
})

setMethod("get_precision", "BenchmarkResults", function(object) {
  precisions <- sapply(1:length(object@pipeline_names), function(i) {
    return(length(object@true_positives[[i]]) / (length(object@true_positives[[i]]) + length(object@false_positives[[i]])))
  })
  names(precisions) <- object@pipeline_names
  return(precisions)
})

setMethod("get_f1", "BenchmarkResults", function(object) {
  recalls <- get_recall(object)
  precisions <- get_precision(object)
  f1s <- sapply(1:length(object@pipeline_names), function(i) {
    return(2 * (precisions[i] * recalls[i]) / (precisions[i] + recalls[i]))
  })
  names(f1s) <- object@pipeline_names
  return(f1s)
})
