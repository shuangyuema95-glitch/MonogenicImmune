#' UMAP Visualization for Batch-Corrected Integrations
#' @title BatchCoordinate_Umap: Plot UMAPs for Integration Methods
#' @description Generate combined UMAP plots colored by sample identity and clustering results
#'   for 5 common single-cell integration methods (Harmony, BBKNN, Scanorama, scVI, SCALEX).
#'   Automatically arranges plots and exports to a single PDF file.
#' @param input Character. Path to the result directory containing CSV coordinate files.
#' @param pcaNumber Numeric. Number of PCA components used for integration (used in file naming).
#' @param fileCoord Character. Name of the CSV file containing UMAP coordinates and metadata.
#' @returnType NULL
#' @return No object returned. A PDF file with UMAP plots is saved to the input directory.
#'   Output includes:
#'   - Paired UMAP plots (colored by samples + clusters) for each integration method
#'   - Combined PDF with organized layout
#' @author MSY

library(ggplot2)
library(patchwork)
library(cowplot)


input="E:\\AID cohort\\code\\result\\"
#fileCoord="AllUMAPLeiden_annMerge.csv"
#pcaNumber=50

my_colors <- c(
  "#CC6666", "#6699CC", "#008280", "#BB0021", "#631879", "#E64B35",
  "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#4A6FE3", "#B5BBE3",
  "#D6BCC0", "#BB7784", "#ED645A", "#8CD17D", "#86C7B4", "#9CD2ED", "#B07AA1",
  "#FF9D9A", "#DAA51B", "#F28E2B", "#D4DE9C", "#8DD593", "#B6992D", "#E58606",
  "#499894", "#FFBE7D", "#FABFD2", "#023FA5", "#CC95C0", "#7F7FD5", "#F3F9A7",
  "#9D7660", "#A0CBE8", "#134E5E", "#4AC29A", "#A8E063", "#F0E442", "#FEC260",
  "#EF9708", "#EAD3C6", "#ADE87C", "#9F8D89", "#DF5734", "#A44E89",
  "#B7DEEA", "#C4DAEC", "#C4612F", "#61BADA", "#7382BC", "#0FCFC0",
  "#8D689D", "#83AB8E", "#D8A0C0", "#ECE399", "#D33F6A", "#B95055",
  "#B85292", "#405993", "#64A776", "#A992C0", "#EA9994", "#007B7F", "#DF75AE", "#007ABA",
  "#A8C7E9", "#00B7CA", "#F37121", "#8CCA7C", "#8B58A4", "#FF5252", "#448AFF",
  "#00ACC1", "#FF9800", "#7B1FA2", "#FF4081", "#1E88E5", "#009688", "#5C6BC0",
  "#FF7043", "#5E35B1", "#43A047", "#FDD835", "#1976D2", "#FFA726",
  "#7E57C2", "#26A69A", "#FFEB3B", "#EC407A", "#3949AB", "#66BB6A", "#FFCA28",
  "#0288D1", "#8E24AA", "#42A5F5", "#AB47BC", "#29B6F6", "#689F38", "#FF5722",
  "#9C27B0", "#1DE9B6", "#994F00", "#8F7700", "#352A86", "#348DB2", "#E06377",
  "#FFB547", "#4E84C4", "#C4AD66", "#77BEDB", "#52BCA3", "#CC99C9",
  "#9ECAE1", "#CCCC99", "#F6B5C8", "#7ED379", "#78C679", "#D65F5F", "#8A9FD1",
  "#C6E090", "#7B685D", "#D9A404", "#998EC3", "#F1A340", "#661E01", "#B1DA99",
  "#A23B72", "#F18D9C", "#593C8C", "#C86070", "#9B5A56", "#727D8C", "#DC9497"
)

BatchCoordinate_Umap<-function(input,pcaNumber,fileCoord){
pcaNumber=pcaNumber
df<-as.data.frame(read.csv(paste0(input,pcaNumber,"_",fileCoord)))
colnames(df)<-gsub("bbknn","BBKNN",colnames(df))
methods <-c("harmony","BBKNN","scanorama","scVI","scalex")

plot_umap_gg <- function(df,methodBatch,
                         pt_size = 0.04) {
  
  df0 <- df[,c("cell_id","samples",grep(methodBatch,colnames(df),value = T))]

  colnames(df0) <-c("cell_id","samples","Cluster","UMAP_1","UMAP_2")
  df0$Cluster<-as.character(df0$Cluster)
  df0$samples<-as.character(df0$samples)
  
  
  p_Cluster <- ggplot(df0, aes(x = UMAP_1, y = UMAP_2, color = Cluster)) +
    geom_point(size = pt_size, shape = 16, stroke = 0) +
    scale_color_manual(values = my_colors) +
    theme_classic() +
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_line(colour = "black",linewidth = 0.3),
      legend.position = "right", aspect.ratio = 1) +
    ggtitle(paste0("Cell count: ", nrow(df0)," | Clusters: ", length(unique(df0$Cluster)))) +
    theme(
      plot.title = element_text(size = 10),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 6)
    ) +guides(color = guide_legend(ncol = 5, override.aes = list(size = 2.5)))
  
  p_samples <- ggplot(df0, aes(x = UMAP_1, y = UMAP_2, color = samples)) +
    geom_point(size = pt_size, shape = 16, stroke = 0) +
    scale_color_manual(values = my_colors) +
    theme_classic() +
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_line(colour = "black",linewidth = 0.3),
      legend.position = "right", aspect.ratio = 1) +
    ggtitle(paste0("Cell count: ", nrow(df0)," | Samples: ", length(unique(df0$samples)))) +
    theme(
      plot.title = element_text(size = 10),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 6)
    ) +guides(color = guide_legend(ncol = 5, override.aes = list(size = 2.5)))
  
 
  return(p_samples + p_Cluster +plot_annotation(
  title = methodBatch,
  theme = theme(plot.title = element_text(size=14, face="bold", hjust=0.5)))+
  plot_layout(ncol = 2))
}

allcoordi<-lapply(methods,function(m){
plot_umap_gg(methodBatch = m,df = df,pt_size = 0.04)})

pdf(file.path(input, "Integration_UMAP_Plots.pdf"), width = 14, height = 10)
for (i in seq(1, length(allcoordi), by = 3)) {
  plots_sub <- allcoordi[i:min(i + 2, length(allcoordi))]
  page <- plot_grid(plotlist = plots_sub, ncol = 1)
  print(page)
}
dev.off()
print("======UMAP plots after removing batch has been done======") 
}

BatchCoordinate_Umap(input=input,pcaNumber = 50,fileCoord = "AllUMAPLeiden_annMerge.csv")
