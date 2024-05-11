#' Summarize (mean, sd) the main stats from a benchmark
#' 
#' @param vc_stats A data frame obtained from the run_vc_benchmark function
#' 
#' @return A summary of the main statistics

summarize_vc_stats <- function(vc_stats) {
  vc_stats %>% 
    group_by(pipeline_names) %>%
    summarize(tp_mean = mean(tp_count),
              tp_sd = sd(tp_count),
              fp_mean = mean(fp_count),
              fp_sd = sd(fp_count),
              fn_mean = mean(fn_count),
              fn_sd = sd(fn_count),
              precision_mean = mean(precision),
              precision_sd = sd(precision),
              recall_mean = mean(recall),
              recall_sd = sd(recall),
              f1_mean = mean(f1),
              f1_sd = sd(f1))
}