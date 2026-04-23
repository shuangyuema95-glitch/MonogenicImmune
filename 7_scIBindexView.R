#' Visualization of scib metrics
#' @title BenchPlots: draw bubble plot and bar plot for scib-metrics
#' @description Perform QC filtering, HVG selection, PCA, and doublet detection for multiple 10X scRNA-seq samples. 
#'   Adds mitochondrial and ribosomal gene percentages, log10GenesPerUMI, and scDblFinder classification to metadata.
#' @param nunmber_var Numeric, the number of methods for removing batch effects
#' @param input Character. Path to parent folder containing csv file of SCIB benchmarking results
#' @returnType ggplot
#' @return A ggplot object contains bubble plot and bar plots
#' @author MSY

library(ggplot2)
library(patchwork)
library(reshape2)
library(dplyr)

BenchPlots<-function(input,number_var){
  
  
  
  Bench<-read.csv(paste0(input,"50_All_Benchmark_Final.csv"))
  rownames(Bench)<-Bench$Embedding;Bench$Embedding<-NULL
  colnames(Bench)<-gsub("\\."," ",colnames(Bench))
  
  
  bio_cols = c("Isolated labels","KMeans NMI","KMeans ARI","Silhouette label","cLISI")
  batch_cols = c("BRAS","iLISI","KBET","Graph connectivity","PCR comparison")
  all_bubble_cols = c(bio_cols, batch_cols)
  
  bubble_color = colorRampPalette(c("#7a2790","#d8b8e3","#e6f4e1","#288b39"))(100)
  bubble_data = Bench[,all_bubble_cols]
  bubble_data$Method = rownames(bubble_data)
  bubble_long = melt(bubble_data, id.vars="Method", variable.name="Metric", value.name="Score")
  Bench$Estimate<-Bench$`Batch correction`+Bench$`Bio conservation`+Bench$Total
  
  bubble_long$Method = factor(bubble_long$Method, levels = rev(Bench%>%arrange(desc(Estimate))%>%rownames()))
  bubble_long$Metric = factor(bubble_long$Metric, levels=all_bubble_cols)
  bubble_long$Score_label <- round(bubble_long$Score,2)
  
  p_bubble <- ggplot(bubble_long,aes(x=Metric,y=Method,fill=Score,size=Score))+
    geom_point(shape=19, aes(color=Score)) + 
    geom_text(aes(label=Score_label),size=3,color="black")+
    scale_color_gradientn(colours=bubble_color,limits=c(0,1))+ 
    scale_fill_gradientn(colours=bubble_color,limits=c(0,1))+   
    scale_size_continuous(range=c(3,11),limits=c(0,1))+
    geom_hline(yintercept=seq(0.5,15,1)[-1][1:(number_var-1)],linetype="dashed",color="gray50")+
    theme_bw() +
    theme(panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA),
          panel.grid = element_blank(),
          axis.text.x=element_text(angle=45,hjust=1,color="black"),
          axis.text.y = element_text(color = "black"),
          legend.position="none",
          axis.ticks = element_line(color = "black"),
          plot.title=element_text(vjust=-1),
          plot.margin=margin(10,5,5,5)) +
    xlab(NULL) + ylab(NULL)
  
  
  bar_color <- colorRampPalette(c("#a2e0c1","#70c8b8","#2774ae","#3949AB"))(100)
  bar_data = Bench[,c("Batch correction","Bio conservation","Total")]
  bar_data$Method = rownames(bar_data)
  bar_long = melt(bar_data, id.vars="Method", variable.name="Metric", value.name="Score")
  bar_long$Metric = factor(bar_long$Metric, levels=c("Batch correction","Bio conservation","Total"))
  bar_long$Method = factor(bar_long$Method, levels = rev(Bench%>%arrange(desc(Estimate))%>%rownames()))
  
  p1 <- ggplot(bar_data, aes(y = Method)) +
    geom_col(aes(x = `Batch correction`, fill = `Batch correction`), width = 0.75) +
    geom_text(aes(x = `Batch correction`, label = round(`Batch correction`,2)), 
              hjust = 1.1, color = "white", size = 3) +
    scale_fill_gradientn(colours = bar_color, limits = c(0,1)) +
    scale_x_continuous(limits = c(0, 1), expand = c(0,0)) +
    theme_classic() + 
    ylab(NULL)+xlab("Batch\ncorrection")+
    theme(legend.position = "none",
          axis.title.x = element_text(angle=45, hjust=1, size=10),
          axis.ticks=element_blank(),
          axis.text.y = element_blank(), axis.line.y = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.line.x = element_blank(), plot.margin = margin(0,2,0,2))
  
  p2 <- ggplot(bar_data, aes(y = Method)) +
    geom_col(aes(x = `Bio conservation`, fill = `Bio conservation`), width = 0.75) +
    geom_text(aes(x = `Bio conservation`, label = round(`Bio conservation`,2)), 
              hjust = 1.1, color = "white", size = 3) +
    scale_fill_gradientn(colours = bar_color, limits = c(0,1)) +
    scale_x_continuous(limits = c(0, 1), expand = c(0,0)) +
    theme_classic() +
    ylab(NULL)+xlab("Bio\nconservation")+
    theme(legend.position = "none",
          axis.title.x = element_text(angle=45, hjust=1, size=10),
          axis.ticks=element_blank(),
          axis.text.y = element_blank(), axis.line.y = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.line.x = element_blank(), plot.margin = margin(0,2,0,2))
  
  p3 <- ggplot(bar_data, aes(y = Method)) +
    geom_col(aes(x = Total, fill = Total), width = 0.75) +
    geom_text(aes(x = Total, label = round(Total,2)), 
              hjust = 1.1, color = "white", size = 3) +
    scale_fill_gradientn(colours = bar_color, limits = c(0,1)) +
    scale_x_continuous(limits = c(0, 1), expand = c(0,0)) +
    theme_classic() +
    ylab(NULL)+xlab("Total")+
    theme(legend.position = "none",
          axis.title.x = element_text(angle=45, hjust=1, size=10),
          axis.ticks=element_blank(),
          axis.text.y = element_blank(), axis.line.y = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.line.x = element_blank(), plot.margin = margin(0,2,0,2))
  
  p_bar <- p1 + p2 + p3 + plot_layout(widths = c(1,1,1))
  
  p_bubble + p_bar + plot_layout(widths = c(1.7, 0.87))
  
}
BenchPlots(input,4)
