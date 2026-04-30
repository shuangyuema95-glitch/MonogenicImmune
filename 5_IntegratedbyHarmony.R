#' Quality Control for multiple scRNA-seq samples
#' @title PrefilterRaw: QC and Doublet Annotation for scRNA-seq
#' @description Perform QC filtering, HVG selection, PCA, and doublet detection for multiple 10X scRNA-seq samples. 
#'   Adds mitochondrial and ribosomal gene percentages, log10GenesPerUMI, and scDblFinder classification to metadata.
#' @param samples Character vector. Sample names corresponding to folders containing 10X matrices.
#' @param fold2 Character. Path to parent folder containing all sample subfolders.
#' @returnType List
#' @return A named list, one element per sample, each containing:
#'   - UMAP plots: for different stages, raw data, scale data 
#'   - Lisi score and boxplots: devtools::install_github("immunogenomics/lisi")
#'   - PBMC1: a combined Seurat by Seurat CAA method
#' @author MSY


library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(ggplot2)
library(lisi)
library(harmony)

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
samples=c("137ACP5","84ADA2","159ADA2","163COPA","88GNAI2","78IL1R1","79IL1R1_2","156ISG15","157LACC1","90OGFRL1","36PLCG1","37PLCG2","108PLD4","109PLD4","110PLD4","155PSMD12","138RIPK1","139RIPK1","140RIPK1","141RIPK1","152RNF31","149SOCS1","162STING","69TBXAS1","70TBXAS1","89TLR7","82TNFAIP3","83TNFAIP3","154TRAF3","151USP18","134ELF4","115KRAS","81NLRC4","91TLR8","148UBA1","1IRAK2","3IRAK2","165TREX1","168TREX1","169TREX1","175IFIH1","176IFIH1","177IFIH1","178ISG15","179STAT4","180SLC7A7","181NOD2","182LPIN2","12HC","17HC","22HC","75HC","76HC","77HC","158HC","160HC","161HC","9HC","HC_GSE139324_GSM4138162","HC_GSE139324_GSM4138163","HC_GSE139324_GSM4138164","HC_GSE139324_GSM4138165","HC_GSE139324_GSM4138166","HC_GSE139324_GSM4138167","HC_10XGenomics10K","HC_GSE156989_GSM4749756","HC_GSE156989_GSM4749762","HC_GSE156989_GSM4749768","HC_GSE199445_GSM5973143","HC_GSE199445_GSM5973144","HC_GSE199445_GSM5973145","HC_GSE199445_GSM5973146","HC_GSE168732_GSM5160432","HC_GSE168732_GSM5160434","HC_GSE168732_GSM5160435","SLE_GSE156989_GSM4749774","SLE_GSE156989_GSM4749779","SLE_GSE156989_GSM4749784","SLE_GSE142016_GSM4217718","SLE_GSE142016_GSM4217719","SLE_GSE142016_GSM4217720","SLE_GSE263931_GSM8207595","SLE_GSE263931_GSM8207597","SLE_GSE263931_GSM8207599","SLE_GSE263931_GSM8207601","SLE_GSE263931_GSM8207605","SLE_GSE263931_GSM8207607","SLE_GSE224198_GSM7017326","SLE_GSE224198_GSM7017329","SLE_GSE224198_GSM7017331","SLE_GSE224198_GSM7017334","JIA_GSE205095_GSM6205132","JIA_GSE205095_GSM6205136","JIA_GSE205095_GSM6205138","JIA_GSE205095_GSM6205140","JIA_GSE205095_GSM6205142","JIA_GSE205095_GSM6205144","KD_GSE168732_GSM5160417","KD_GSE168732_GSM5160420","KD_GSE168732_GSM5160422","KD_GSE168732_GSM5160424","KD_GSE168732_GSM5160427","KD_GSE168732_GSM5160430","STING_GSE226598_GSM7079986","STING_GSE226598_GSM7079989","STING_GSE226598_GSM7079991","PLCG2_Figshare","OTULIN_GSE199445_GSM5973147","OTULIN_GSE199445_GSM5973148","OTULIN_GSE199445_GSM5973149")
input="/public/home/yuxiaomingroup/msy/work/AID/juzhen/clean/result/"
options(future.globals.maxSize = 100 * 1024^3)

Seu2Merge_Lisi<-function(input,samples,listtype="top3000",pcaNumber=NULL){

  pcaNumber=pcaNumber
  targetfold=paste0(input,pcaNumber)
  if (!dir.exists(targetfold)) {dir.create(targetfold, recursive = TRUE)}
  
  plot_umap_gg <- function(Seurat_obj,reduction = "",group_by = "samples",
                         pt_size = 0.04,sufix = "",pname="",pcaNumber=pcaNumber) {
  
  
  df <- as.data.frame(Seurat_obj@reductions[[reduction]]@cell.embeddings)
  colnames(df) <- c("UMAP_1", "UMAP_2")
  df$Cluster<-Seurat_obj@meta.data$samples
  
  p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = Cluster)) +
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
    ggtitle(paste0("Cell count: ", ncol(Seurat_obj)," | Samples: ", length(unique(df$Cluster)))) +
    theme(
      plot.title = element_text(size = 10),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 6)
    ) +guides(color = guide_legend(ncol = 5, override.aes = list(size = 2.5)))
  
  filename <- file.path(targetfold, paste0(pcaNumber, sufix, "_", pname, "UMAP.pdf"))
  ggsave(filename, plot = p, width = 12.6, height = 11.5)}
 
  
  pbmc1 <-LoadH5Seurat(paste0(input,"pbmcmerged.h5Seurat"))
  pbmc1 <- NormalizeData(pbmc1, normalization.method = "LogNormalize", scale.factor = 10000)

  if(listtype=="top3000"){
    cc_genes <- as.character(c(cc.genes$s.genes,cc.genes$g2m.genes))
    hvgs_top3000<-as.character(readRDS(paste0(input,"hvgs_top3000.rds")))
    integ_features <- hvgs_top3000[!hvgs_top3000%in%cc_genes]
    VariableFeatures(pbmc1) <-integ_features
    sufix<-"t3k"
  }else{
    cc_genes <- as.character(c(cc.genes$s.genes,cc.genes$g2m.genes))
    mergedhvg<-as.character(readRDS(paste0(input,"merge_hvgs.rds")))
    integ_features <- mergedhvg[!mergedhvg%in%cc_genes]
    
    VariableFeatures(pbmc1) <-integ_features
    sufix<-"mhvg"
  }

  ## (1) raw
  pbmc2 <- ScaleData(pbmc1)
  pbmc2 <- RunPCA(pbmc2,npcs = pcaNumber)
  pbmc2 <- RunUMAP(pbmc2,reduction = "pca",dims = 1:pcaNumber)
  
  cells_use <- sample(Cells(pbmc2), size = 30000)
  pbmc2_sub <- subset(pbmc2, cells = cells_use)
  lisi_A<-compute_lisi(Embeddings(pbmc2_sub,'pca'),pbmc2_sub@meta.data,'samples')
  rm(cells_use);rm(pbmc2_sub);gc()
  
  plot_umap_gg(pbmc2,reduction ='umap',sufix = sufix,pname="ggRaw",pcaNumber=pcaNumber)
  rm(pbmc2);gc()
  print("---------raw UMAP has been done and saved to pdf------------\n")
    
  ## (2) Scale
  g2m_genes <- CaseMatch(search = cc.genes$g2m.genes, match = rownames(pbmc1))
  s_genes  <- CaseMatch(search = cc.genes$s.genes,  match = rownames(pbmc1))
  pbmc1 <- CellCycleScoring(pbmc1, g2m.features = g2m_genes, s.features = s_genes)
  pbmc1 <- ScaleData(pbmc1, vars.to.regress = c("percent_mt", "S.Score", "G2M.Score",'nCount_RNA'))
  pbmc1 <- RunPCA(pbmc1, features = VariableFeatures(pbmc1), npcs = pcaNumber)
  pbmc1 <- RunUMAP(pbmc1,reduction = "pca",dims = 1:pcaNumber)
  
  cells_use <- sample(Cells(pbmc1), size = 30000)
  pbmc1_sub <- subset(pbmc1, cells = cells_use)
  lisi_B<-compute_lisi(Embeddings(pbmc1_sub,'pca'),pbmc1_sub@meta.data,'samples')
  rm(cells_use);rm(pbmc1_sub);gc()
  
  plot_umap_gg(pbmc1,reduction ='umap',sufix = sufix,pname="ggScaled",pcaNumber=pcaNumber)
  print("---------scaled2 UMAP has been done and saved to pdf------------\n")

  ## (3) Run Harmony using pbmc1, since pbmc2 was run twice with scaling.
  gc()
  start_total <- proc.time()
  before_mem <- sum(gc()[, 1])
  pbmc1 <- RunHarmony(pbmc1, group.by.vars = "samples", dims.use = 1:pcaNumber, theta = 2) 
  end_total <- proc.time()
  after_mem <- sum(gc()[, 1])
  Runtime_Mins <- round((end_total - start_total)[["elapsed"]] / 60, 3)
  CPU_Time_Seconds <- round((end_total - start_total)[["user.self"]], 2)  
  Memory_Used_GB <- round((after_mem - before_mem) / 1024, 3)
  
  harmony_resource <- data.frame(Step = "Harmony",Runtime_Mins,CPU_Time_Seconds, Memory_Used_GB)
  openxlsx::write.xlsx(harmony_resource,file=file.path(targetfold, paste0(sufix,pcaNumber, "_HarmonyUsage.xlsx")))

  
  pbmc1 <- FindNeighbors(pbmc1, dims = 1:pcaNumber, features =VariableFeatures(pbmc1),reduction = "harmony")
  pbmc1 <- FindClusters(pbmc1, resolution = 0.3)
  set.seed(1)
  pbmc1 <- RunUMAP(pbmc1, reduction = "harmony", dims = 1:pcaNumber, reduction.name = "harmony_umap")
  pbmc1$harmony_clusters <- pbmc1$seurat_clusters
    
  cells_use <- sample(Cells(pbmc1), size = 30000)
  pbmc1_sub <- subset(pbmc1, cells = cells_use)
  lisi_C<-compute_lisi(Embeddings(pbmc1_sub,'harmony'),pbmc1_sub@meta.data,'samples')
  rm(cells_use);rm(pbmc1_sub);gc()
 
  plot_umap_gg(pbmc1,reduction ='harmony_umap',sufix = sufix,pname="ggHarmony",pcaNumber=pcaNumber)
  print("---------Harmony UMAP has been done and saved to pdf------------\n")
  
  
  pbmc1 <- DietSeurat(pbmc1,counts = TRUE,data = TRUE,scale.data = TRUE)
  SaveH5Seurat(pbmc1, filename = file.path(targetfold, paste0(pcaNumber, sufix, "_harmony.h5Seurat")))
  rm(pbmc1);gc()
 
  ## (4) boxplot to compare lisi score across raw, scale and harmony stages
  lisi_data=data.frame(
      step=c(rep("Raw",nrow(lisi_A)),rep("Z-score",nrow(lisi_B)),rep("Harmony",nrow(lisi_C))),
      LISI=c(lisi_A$samples,lisi_B$samples,lisi_C$samples))
  lisi_data$step<-factor(lisi_data$step,levels = c("Raw","Z-score","Harmony"))
  pici_lisi<-ggplot(lisi_data, aes(x = step, y = LISI, fill = step)) +
      geom_boxplot(outlier.size = 0.01, width =0.7, position = position_dodge(width = 0.7)) +  
      stat_summary(aes(group = step, color = step), fun = median, geom = "line", size = 1.2) +  
      scale_fill_manual(values = c("#727D8C", "#023FA5","#3949AB")) +  
      labs(x = "", y = "LISI score", fill = "Step") +
      theme_minimal() +
      theme(axis.text = element_text(color = "black"),
            legend.position = "none",
            axis.ticks.y = element_line(color="black"),
            axis.title.x = element_blank(), 
            panel.spacing = unit(1, "lines"),  
            panel.border = element_rect(color = "black", fill = NA, size = 0.35))
  
  saveRDS(lisi_data, file = file.path(targetfold, paste0(pcaNumber, sufix, "_lisidata.RDS")))
  ggsave(plot = pici_lisi,filename = file.path(targetfold, paste0(pcaNumber, sufix, "_Lisibox.pdf")),width = 6.9, height = 7.5)
  
  rm(list=ls())
  gc()
}

before2time <- Sys.time()
Seu2Merge_Lisi(input,samples,"top3000",pcaNumber=50)
after2time <- Sys.time()
total_run_time <- difftime(after2time, before2time, units = "mins")
cat(total_run_time )
