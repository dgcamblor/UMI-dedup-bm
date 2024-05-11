#' Generate a comparative barplot for deduplication benchmarking based on the provided summarized data.
#' 
#' @param summarized_data The summarized data used for generating the barplot.
#' @param stat It can be either "tp" for true positives or "fp" for false positives.
#' @param level_order The order of the pipeline names in the barplot. 
#' @param full_names The full names corresponding to the pipeline names. 
#' @param title_lab The title of the barplot.
#' @param image_path The file path where the generated barplot image will be saved.
#' @param hjust The horizontal justification for the text labels on the bars.
#' @param vjust The vertical justification for the text labels on the bars.
#' 
#' @return The generated barplot as a ggplot object.

dedup_bm_barplot <- function(summarized_data, stat, level_order, full_names, title, image_path, hjust = -0.55, vjust = 0.5) {
  if (length(level_order) != length(pipeline_names)) {
    print("Check the pipeline names and the level order")
    stop()
  }

  aes_y <- ifelse(stat == "tp", "tp_mean", "fp_mean")
  aes_sd <- ifelse(stat == "tp", "tp_sd", "fp_sd")
  y_lab <- ifelse(stat == "tp", "True positives", "False positives")
  
  data$pipeline_names <- factor(data$pipeline_names, levels = level_order, labels = full_names) 
  
  tiff(image_path, units = "px", width = 2828, height = 2044, res = 400)
  plot <- ggplot(data, aes(x = pipeline_names, y = !!sym(aes_y), fill = pipeline_names)) +
    geom_bar(stat = "identity", color = "black") +
    geom_errorbar(aes(ymin = !!sym(aes_y) - !!sym(aes_sd), ymax = !!sym(aes_y) + !!sym(aes_sd)), width = 0.2) +
    labs(x = "", y = y_lab, fill = "Deduplication", title = title) +
    geom_text(aes(label = sprintf("%.2f", !!sym(aes_y)), y = 0), hjust = hjust, vjust = vjust, position = position_dodge(width = 0.9), size = 6, angle = 90) +
    theme_bw() +
    theme(title = element_text(face = "bold"))
  print(plot)
  dev.off()
  
  return(plot)
}