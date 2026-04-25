#' LISI Integration Evaluation for scRNA-seq Integration Workflow
#' @title LISI_compare: PCA Number Optimization & Integration Performance Evaluation
#' @description Comprehensive evaluation of scRNA-seq integration performance across different PCA dimensions (30/35/40/45/50).
#'   Automatically loads LISI score results for two HVG selection strategies (top3k HVG & merged HVG), generates standardized
#'   facet boxplots for LISI scores, calculates integration improvement metrics, and creates square heatmaps for visualization.
#'   Used to determine optimal PCA number for downstream integration analysis.
#' @param input Character. Path to the working directory containing PCA-numbered LISI result RDS files.
#' @returnType No direct return (automatically saves plots to input directory).
#' @return Two PDF figures saved to input directory:
#'   - hvg_mer_boxlisi.pdf: Faceted boxplots of LISI scores across PCA numbers & integration steps (Raw/Z-score/Harmony)
#'   - hvg_mer_diffSum.pdf: Square heatmaps of mean LISI values, differences, and total integration scores
#' @author MSY


library(tidyverse)
library(patchwork)
input="E:\\AID cohort\\code\\result\\PCAnumber"

LISI_compare<-function(input){
setwd("input")
hvgt3_lisi<-data.frame()
for(i in c(30,35,40,45,50)){
  s<-readRDS(paste0(i,"t3k_lisidata.RDS"))
  print(s%>%group_by(step)%>%summarise(junzhi=median(LISI)))
  s$PCAnumber<-i
  hvgt3_lisi<-rbind(hvgt3_lisi,s)
}

hvgmer_lisi<-data.frame()
for(i in c(30,35,40,45,50)){
  s<-readRDS(paste0(i,"mhvg_lisidata.RDS"))
  print(s%>%group_by(step)%>%summarise(junzhi=mean(LISI)))
  s$PCAnumber<-i 
  hvgmer_lisi<-rbind(hvgmer_lisi,s)
}



Box_View<-function(data1){
  return(ggplot(data1, aes(x=step, y=LISI, fill=step)) +
           geom_boxplot(position=position_dodge(0.9), width=0.8, outlier.size=0.01, key_glyph="rect") +
           stat_summary(aes(group=step), fun=median, geom="line", size=1.2, color="black") +
           scale_fill_manual(values=c("#727D8C", "#B5BBE3","#1976D2")) +
           facet_wrap(~paste0("PCA: ", PCAnumber), nrow=1) +
           labs(x="", y="LISI score", fill="Step") +
           theme_minimal() +
           theme(
             axis.text=element_text(color="black", size=11),
             axis.text.x=element_blank(),
             axis.ticks.x=element_blank(),
             axis.ticks.y=element_line(color="black"),
             legend.position="right",
             legend.key.size=unit(0.32, "cm"), 
             legend.text=element_text(size=8), 
             legend.title=element_text(size=9),
             panel.border=element_rect(color="black", fill=NA, size=0.35),
             strip.background=element_rect(fill="gray90", color="black"),
             strip.text=element_text(size=10),
             plot.margin=margin(10,10,10,10)))
}
P1<-Box_View(hvgt3_lisi)
P2<-Box_View(hvgmer_lisi)
Diff_Cal<-function(data1){

lisi_mean <- data1 %>%
  group_by(step, PCAnumber) %>%
  summarise(meanValue = mean(LISI), .groups = "drop")

lisi_wide <- lisi_mean %>%
  pivot_wider(names_from = step, values_from = meanValue)

lisi_result <- lisi_wide %>%
  mutate(
    diff_Harmony_Raw = Harmony - Raw,
    diff_Harmony_Zscore = Harmony - `Z-score`,
    total_score = diff_Harmony_Raw + diff_Harmony_Zscore
  ) %>%arrange(desc(total_score)) %>%as.data.frame()

colnames(lisi_result)<-c(
  "PCA_number","raw","Z-score",
  "Harmony","Harmony-raw","Harmony-Zscore","total diffence")

lisi_result2<-lisi_result %>% 
  mutate(PCA_number = factor(PCA_number, levels = rev(lisi_result$`PCA_number`))) %>% 
  pivot_longer(cols = -PCA_number, names_to = "metric", values_to = "value")

lisi_result2$metric<-factor(
lisi_result2$metric,levels = c("raw","Z-score","Harmony",
"Harmony-raw","Harmony-Zscore","total diffence"))
return(ggplot(
  lisi_result2,
  aes(x = metric, y = PCA_number)) +
  geom_tile(aes(fill = value), color = "black", size = 0.4) +
  geom_text(aes(label = round(value, 2)), color = "black", size = 3.6) +
  scale_fill_gradient(low = "#B5BBE3", high = "#1976D2") +
  labs(x = "", y = "PCA Number", fill = "Value") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(color = "black", size = 8),#angle = 45, hjust = 1, 
    axis.text.y = element_text(color = "black", size = 10),
    panel.grid = element_blank(),
    legend.position = "right"
))



}
P3<-Diff_Cal(hvgt3_lisi)
P4<-Diff_Cal(hvgmer_lisi)
boxplots<-P1/P2
diffBubble<-(P3+theme(legend.position = "none"))/(P4+theme(legend.position = "none"))
ggsave(plot = boxplots,filename = file.path(input,  "hvg_mer_boxlisi.pdf"),width = 8.9, height = 4.5)
ggsave(plot = diffBubble,filename = file.path(input,  "hvg_mer_diffSum.pdf"),width = 5.9, height = 4.5)

}

