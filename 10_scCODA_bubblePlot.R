#' scCODA Bubble Plot for Cell Type Effect Visualization
#' @title plot_scCODA_bubble: Visualize scCODA Results as a Bubble Plot
#' @description Generate a publication-ready bubble plot from scCODA differential cell 
#'   composition analysis. Plot shows log2-fold changes colored by effect size, 
#'   with point size representing statistical significance.
#' @param csv_file Character. Path to the CSV file containing scCODA results.
#' @return A ggplot2 object displaying the bubble plot.
#' @author MSY
#' @export

library(ggplot2)
library(readr)

plot_scCODA_bubble <- function(csv_file) {
  df <- read_csv(csv_file, show_col_types = FALSE)
  ggplot(df, aes(x = Group, y = `Cell Type`)) +
    geom_point(aes(color = log2FC, size = Significant)) +
    scale_color_gradient2(
      low = "#4A148C", mid = "#FFFFFF", high = "#E64A19",
      limits = c(-2, 1.5),
      name = "Effect\nlog2(FC)"
    ) +
    scale_size_manual(
      values = c("TRUE" = 6, "FALSE" = 2.5),
      name = "Significance Effect"
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(color = "black", size = 10),
      axis.ticks = element_line(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      legend.title = element_text(color = "black", size = 10),
      legend.text = element_text(color = "black")
    ) +
    labs(x = "", y = "")
}

plot_scCODA_bubble("scCODA_all_results_bubble.csv")