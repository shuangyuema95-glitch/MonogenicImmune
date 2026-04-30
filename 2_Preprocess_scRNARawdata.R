#' Quality Control for multiple scRNA-seq samples
#' @title PrefilterRaw: QC and Doublet Annotation for scRNA-seq
#' @description Perform QC filtering, HVG selection, PCA, and doublet detection for multiple 10X scRNA-seq samples. 
#'   Adds mitochondrial and ribosomal gene percentages, log10GenesPerUMI, and scDblFinder classification to metadata.
#' @param samples Character vector. Sample names corresponding to folders containing 10X matrices.
#' @param fold2 Character. Path to parent folder containing all sample subfolders.
#' @param outdir Character. Output directory for QC results (default: current working directory).
#' @returnType List
#' @return A named list, one element per sample, each containing:
#'   - QC_seurat: Seurat object after QC, normalization, HVG selection, PCA, and scDblFinder annotation
#'   - summary: data.frame summarizing cell and gene counts before and after QC, median nFeature_RNA, nCount_RNA, percent.mt, and percent.ribo
#' @author MSY

setwd("E:\\AID cohort\\code")

library(Seurat)
library(data.table)
library(tidyverse)
library(scDblFinder)
library(SingleCellExperiment)
library(SeuratDisk)
library(openxlsx)


metainfo<-as.data.frame(fread("E:\\AID cohort\\code\\metainfo.txt",header=TRUE))
QCdata<-data.frame(
  metric=c("3' V3","3' V4","3' V2","5' V1","5' V2","Drop-seq"),
  mt_threshold=c(25,30,20,20,20,20),
  rb_threshold=c(60,60,65,60,60,60),
  feature_threshold=c(6000,8000,4000,6000,6000,4000),
  count_threshold=c(50000,60000,25000,50000,40000,25000)
)

samples=c("137ACP5","84ADA2","159ADA2","163COPA","88GNAI2","78IL1R1","79IL1R1_2","156ISG15","157LACC1","90OGFRL1","36PLCG1","37PLCG2","108PLD4","109PLD4","110PLD4","155PSMD12","138RIPK1","139RIPK1","140RIPK1","141RIPK1","152RNF31","149SOCS1","162STING","69TBXAS1","70TBXAS1","89TLR7","82TNFAIP3","83TNFAIP3","154TRAF3","151USP18","134ELF4","115KRAS","81NLRC4","91TLR8","148UBA1","1IRAK2","3IRAK2","165TREX1","168TREX1","169TREX1","175IFIH1","176IFIH1","177IFIH1","178ISG15","179STAT4","180SLC7A7","181NOD2","182LPIN2","12HC","17HC","22HC","75HC","76HC","77HC","158HC","160HC","161HC","9HC","HC_GSE139324_GSM4138162","HC_GSE139324_GSM4138163","HC_GSE139324_GSM4138164","HC_GSE139324_GSM4138165","HC_GSE139324_GSM4138166","HC_GSE139324_GSM4138167","HC_10XGenomics10K","HC_GSE156989_GSM4749756","HC_GSE156989_GSM4749762","HC_GSE156989_GSM4749768","HC_GSE199445_GSM5973143","HC_GSE199445_GSM5973144","HC_GSE199445_GSM5973145","HC_GSE199445_GSM5973146","HC_GSE168732_GSM5160432","HC_GSE168732_GSM5160434","HC_GSE168732_GSM5160435","SLE_GSE156989_GSM4749774","SLE_GSE156989_GSM4749779","SLE_GSE156989_GSM4749784","SLE_GSE142016_GSM4217718","SLE_GSE142016_GSM4217719","SLE_GSE142016_GSM4217720","SLE_GSE263931_GSM8207595","SLE_GSE263931_GSM8207597","SLE_GSE263931_GSM8207599","SLE_GSE263931_GSM8207601","SLE_GSE263931_GSM8207605","SLE_GSE263931_GSM8207607","SLE_GSE224198_GSM7017326","SLE_GSE224198_GSM7017329","SLE_GSE224198_GSM7017331","SLE_GSE224198_GSM7017334","JIA_GSE205095_GSM6205132","JIA_GSE205095_GSM6205136","JIA_GSE205095_GSM6205138","JIA_GSE205095_GSM6205140","JIA_GSE205095_GSM6205142","JIA_GSE205095_GSM6205144","KD_GSE168732_GSM5160417","KD_GSE168732_GSM5160420","KD_GSE168732_GSM5160422","KD_GSE168732_GSM5160424","KD_GSE168732_GSM5160427","KD_GSE168732_GSM5160430","STING_GSE226598_GSM7079986","STING_GSE226598_GSM7079989","STING_GSE226598_GSM7079991","PLCG2_Figshare","OTULIN_GSE199445_GSM5973147","OTULIN_GSE199445_GSM5973148","OTULIN_GSE199445_GSM5973149")
length(samples)
fold2="E:\\AID cohort\\code"
outdir="E:\\AID cohort\\code\\result"


PrefilterRaw<-function(samples,fold2,outdir){

  message("=====QC check start====")
  start_time <- Sys.time()  
  if (!dir.exists(outdir)) dir.create(outdir)
  
  # metric	mt	Feature	Count
  # 3‘ V3	<25%	<6000	<50000
  # 3‘ V4	<30%	<8000	<60000
  # 3‘ V2	<20%	<4000	<25000
  # 5’ V1	<20%	<6000	<50000
  # 5’ V2	<20%	<6000	<40000
  # Drop-seq	<20%	<4000	<25000
  
 
  res_list<-lapply(samples,function(x){
  message("\n--- Processing sample: ", x, " ---")
  
  ### genes detected in fewer than 3 cells were removed  
  counts <- Read10X(paste(fold2, x, sep="/"))
  beforeQC_n <- ncol(counts)  
  beforeQC_g <- nrow(counts)   
  s <- CreateSeuratObject(counts = counts, min.cells = 3)

  s$samples<-x
  s[["percent.mt"]] <- PercentageFeatureSet(s, pattern = "^MT-")
  s[["percent.ribo"]] <- PercentageFeatureSet(s, pattern = "^RP[SL]")
  s$log10GenesPerUMI <- log10(s$nFeature_RNA) / log10(s$nCount_RNA)
  
  platelet_genes<-c("PF4","PPBP","ITGA2B","GNG11","TUBB1")
  s[["percent.platelet"]] <- PercentageFeatureSet(s, features = intersect(platelet_genes, rownames(s@assays$RNA@counts)))
  hemo_genes<-c("HBB", "HBD", "HBA1", "HBA2")
  s[["percent.rbc"]] <- PercentageFeatureSet(s, features = intersect(hemo_genes, rownames(s@assays$RNA@counts)) )
  
  
  # filter low-quality cells
  thres<-QCdata[QCdata$metric==filter(metainfo,dataset==x)%>%select(chemistry)%>%
       as.matrix()%>%as.character(),2:5]%>%as.matrix()%>%as.numeric()
  s <- subset(s,
                subset = nCount_RNA < thres[4] &
                nFeature_RNA > 200 & nFeature_RNA < thres[3] &
                percent.mt < thres[1] &
                log10GenesPerUMI > 0.7)

  # before find variable genes, some type of counfounding genes were removed
  all.genes <- rownames(s)
  exclude.genes <- unique(c(
    grep("^IG[HKL]", all.genes, value = TRUE),   # BCR（更精确）
    grep("^TRAV", all.genes, value = TRUE),
    grep("^TRAJ", all.genes, value = TRUE),
    grep("^TRAC", all.genes, value = TRUE),
    
    grep("^TRBV", all.genes, value = TRUE),
    grep("^TRBJ", all.genes, value = TRUE),
    grep("^TRBC", all.genes, value = TRUE),
    
    grep("^TRGV", all.genes, value = TRUE),
    grep("^TRGJ", all.genes, value = TRUE),
    grep("^TRGC", all.genes, value = TRUE),
    
    grep("^TRDV", all.genes, value = TRUE),
    grep("^TRDJ", all.genes, value = TRUE),
    grep("^TRDC", all.genes, value = TRUE),
    
    grep("^RPS", all.genes, value = TRUE),
    grep("^RPL", all.genes, value = TRUE),
    grep("^HSP", all.genes, value = TRUE),
    grep("^MT-", all.genes, value = TRUE)
  ))
  keep.genes <- setdiff(all.genes, exclude.genes)
  
  # Normalization and FindVariable for each datasets and before removing doublets
  s <- subset(s, features = keep.genes)
  s <- NormalizeData(s, normalization.method = "LogNormalize", scale.factor = 10000)
  s <- FindVariableFeatures(s, selection.method = "vst", nfeatures = 4000)

  #s <- ScaleData(s, features = hvg_genes)#, vars.to.regress = c("percent.mt", "percent.ribo", "S.Score", "G2M.Score"))
  #s <- RunPCA(s, features = hvg_genes) after merge, still need scale
  
  # Doublet detection
  s1 <- as.SingleCellExperiment(s) # sce_sce <- scDblFinder(sce_sce,dbr=0.03)
  s1<-scDblFinder(s1,dbr=0.03)#identical(colnames(s1),rownames(s@meta.data))
  #singlets <- colnames(s1)[s1$scDblFinder.class == "singlet"]
  s@meta.data$scDblFinder.class<-s1$scDblFinder.class
  #s <- subset(s, cells = singlets);
  rm(s1);gc()
  
  afterQC_n<-length(colnames(s)[s$scDblFinder.class == "singlet"])
  afterQC_g<-nrow(s)
  
  df <- data.frame(
    sample = x,
    BeforeControl_cells = beforeQC_n,
    BeforeControl_genes =beforeQC_g,
    AfterControl_cells = afterQC_n,
    AfterControl_genes =afterQC_g,
    median_nFeature = median(s@meta.data$nFeature_RNA, na.rm = TRUE),
    median_nCount = median(s@meta.data$nCount_RNA, na.rm = TRUE),
    median_mt = median(s@meta.data$percent.mt, na.rm = TRUE),
    median_ribo = median(s@meta.data$percent.ribo, na.rm = TRUE)
  )
  
  setwd(outdir)
  SaveH5Seurat(DietSeurat(s, scale.data = FALSE),
               filename = paste0(x,".h5Seurat"))
  
  
  
  return(list(summary = df))
  
  
  })
  do.call("rbind",res_list)
  
  summary_df <- do.call(rbind, lapply(res_list, function(x) x$summary))
  summary_df$chemistry<-as.character(as.matrix(metainfo[match(summary_df$sample,metainfo$dataset),'chemistry']))
  summary_df<- as.data.frame(lapply(summary_df, as.character), stringsAsFactors = FALSE)
  if (!is.null(outdir)) {
    write.xlsx(summary_df, file = file.path(outdir, "scRNA_AfterQC.xlsx"))
  }
  
  message("=====QC done====")
  end_time <- Sys.time()  
  message("Total runtime: ", round(difftime(end_time, start_time, units = "mins"), 2), " minutes")
}
QCallsamples<-PrefilterRaw(samples,fold2,outdir)






