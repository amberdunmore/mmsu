#' @import ggplot2
#' @importFrom stats time as.formula setNames
#' @importFrom ggplot2 aes facet_wrap
utils::globalVariables(c("value", "variable"))

#' Plot Malaria Model Results
#'
#' This function plots the results of a malaria model simulation.
#'
#' @param model_results A data frame containing the model results or an odin model object.
#' @param output_vars A character vector specifying the output variables to plot. If NULL, all output variables will be plotted.
#'
#' @return A ggplot object.
#' @export
plot_model <- function(model_results, output_vars = NULL) {

  # grab just the result we want
  results_df <- data.frame(time = model_results[, "t"], model_results[, output_vars, drop = FALSE])
  results_long <- tidyr::pivot_longer(results_df, cols = -time, names_to = "variable", values_to = "value")

  # create simple plot
  p <- ggplot() +
    geom_line(data = results_long, aes(x = time, y = value, color = variable), na.rm = TRUE)
  p <- p + labs(x = "Time", y = "Value", color = "Variable")
  p <- p + theme_minimal()

  print(p)
  invisible(p)
}
