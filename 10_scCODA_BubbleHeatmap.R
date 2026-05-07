#' scCODA Bubble Plot for Cell Type Effect Visualization
#' @title plot_scCODA_bubble: Visualize scCODA Results as a Bubble Plot
#' @description Generate a publication-ready bubble plot from scCODA differential cell 
#'   composition analysis. Plot shows log2-fold changes colored by effect size, 
#'   with point size representing statistical significance.
#' @param input Character. Path to the input directory containing the result file.
#' @param csv_file Character. Name of the CSV file with scCODA differential results.
#' @return A ggplot2 object displaying the bubble plot.

library(ggplot2)
library(readr)

input<-"C:\\Users\\Administrator\\Desktop\\scCODA\\"

library(tidyverse)
library(cowplot)

plot_scCODA_bubble <- function(input, csv_file) {
  df <- read_csv(paste0(input, csv_file), show_col_types = FALSE)
  ggplot(df, aes(x = Group, y = `Cell Type`)) +
    geom_point(aes(color = log2FC, size = Significant)) +
    scale_color_gradient2(low = "#4A148C", mid = "#FFFFFF", high = "#E64A19", limits = c(-2, 1.5), name = "log2(FC)") +
    scale_size_manual(values = c("TRUE" = 6, "FALSE" = 2.5), name = "Significant") +
    theme_bw() +
    theme(legend.key.height=unit(0.4,"cm"),
          legend.key.width=unit(0.35,"cm"),
          axis.text = element_text(size = 10, color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank(),
          legend.title = element_text(size = 10),
          legend.text = element_text(color = "black")) +
    labs(x = "", y = "")
}


#' scCODA Heatmap for Significant Change Frequency
#' @title plot_scCODA_heatmap: Visualize scCODA Significant Count Matrix
#' @description Generate a publication-ready heatmap showing the frequency of 
#'   significant cell composition changes across iterative reference scCODA runs.
#' @param input Character. Path to the input directory containing the count matrix.
#' @param csv_file Character. Name of the CSV file with scCODA significant count matrix.
#' @return A ggplot2 object displaying the frequency heatmap.
plot_scCODA_heatmap <- function(input, csv_file) {
  mat <- read.csv(paste0(input, csv_file), row.names = 1, check.names = F)
  mat_long <- mat %>% rownames_to_column("CellType") %>% pivot_longer(-CellType, names_to = "Group", values_to = "Count")
  ggplot(mat_long, aes(x = Group, y = CellType)) +
    geom_tile(aes(fill = Count), color = "white", linewidth = 0.2) +
    scale_fill_gradientn(colors = c("#FFF5F0","#FDBB84","#FC8D59","#E34A33","#B30000","#67000D")) +
    theme_bw() +
    theme(panel.border=element_blank(),
          legend.key.height=unit(0.4,"cm"),
          legend.key.width=unit(0.35,"cm"),
          axis.text = element_text(size = 10, color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank(),
          legend.title = element_text(size = 10),
          legend.text = element_text(color = "black")) +
    labs(x = "", y = "", fill = "Count")
}

# Run and plot
p1 <- plot_scCODA_heatmap("./", "scCODA_heatmap_count_matrix.csv")
p2 <- plot_scCODA_bubble("./", "scCODA_HMC_for_R.csv")
plot_grid(p1, p2, ncol = 2, align = "hv")
