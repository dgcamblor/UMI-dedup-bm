#' Generate a recall-precision plot
#' 
#' @param data Stats from a benchmarking
#' @param pipeline_names Names of the pipelines to be analyzed
#' @param level_order The order of the pipeline names in the barplot. 
#' @param full_names The full names corresponding to the pipeline names. 
#' @param title Title of the plot
#' @param image_path Path to save the image (.tiff)

generate_rp_plot <- function(data, pipeline_names, level_order, full_names, title, image_path) {
  if (length(level_order) != length(pipeline_names)) {
    print("Check the pipeline names and the level order")
    stop()
  }

  # Change pipeline factor
  data$pipeline_names <- factor(data$pipeline_names, levels = level_order, labels = full_names) 
  
  tiff(image_path, units = "px", width = 2828, height = 2044, res = 400)
  plot <- ggplot(data, aes(x = recall, y = precision, color = pipeline_names)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = "Recall", y = "Precision", color = "Deduplication", title = title) +
    theme_bw() +
    theme(title = element_text(face = "bold"))
  print(plot)
  dev.off()
  
  return(plot)
}