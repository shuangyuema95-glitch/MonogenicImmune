#' Quality Control for scRNA-seq samples
#' @title QC for a single scRNA-seq sample
#' @description Generate QC violin plots and summary statistics for a single 10X sample.
#' @param x Character. Sample name or folder containing the 10X matrix.
#' @param fold2 Character. Path to parent folder containing all samples.
#' @param outdir Character. Output directory (default: "QC_output").
#' @returnType List
#' @return A list containing:
#'   - plot: ggplot object of violin plots (nFeature_RNA, nCount_RNA, percent.mt, percent.ribo)
#'   - summary: data.frame of QC metrics (median, mean, 1% and 99% quantiles, cell number)
#' @author MSY


library(Seurat)
library(ggplot2)
library(patchwork)
library(future)
library(future.apply)
library(openxlsx)
library(cowplot)
library(data.table)


setwd("E:\\AID cohort\\code")
samples=c("137ACP5","84ADA2","159ADA2","163COPA","88GNAI2","180SLC7A7","78IL1R1","79IL1R1_2","156ISG15","157LACC1","90OGFRL1","36PLCG1","37PLCG2","108PLD4","109PLD4","110PLD4","155PSMD12","138RIPK1","139RIPK1","140RIPK1","141RIPK1","152RNF31","149SOCS1","162STING","69TBXAS1","70TBXAS1","89TLR7","82TNFAIP3","83TNFAIP3","154TRAF3","151USP18","134ELF4","115KRAS","81NLRC4","91TLR8","148UBA1","1IRAK2","3IRAK2","165TREX1","168TREX1","169TREX1","175IFIH1","176IFIH1","177IFIH1","178ISG15","179STAT4","12HC","17HC","22HC","75HC","76HC","77HC","158HC","160HC","161HC","9HC","HC_GSE139324_GSM4138162","HC_GSE139324_GSM4138163","HC_GSE139324_GSM4138164","HC_GSE139324_GSM4138165","HC_GSE139324_GSM4138166","HC_GSE139324_GSM4138167","HC_10XGenomics10K","HC_GSE156989_GSM4749756","HC_GSE156989_GSM4749762","HC_GSE156989_GSM4749768","HC_GSE199445_GSM5973143","HC_GSE199445_GSM5973144","HC_GSE199445_GSM5973145","HC_GSE199445_GSM5973146","HC_GSE168732_GSM5160432","HC_GSE168732_GSM5160434","HC_GSE168732_GSM5160435","HC_CELLxGENE","SLE_GSE156989_GSM4749774","SLE_GSE156989_GSM4749779","SLE_GSE156989_GSM4749784","SLE_GSE142016_GSM4217718","SLE_GSE142016_GSM4217719","SLE_GSE142016_GSM4217720","SLE_GSE263931_GSM8207595","SLE_GSE263931_GSM8207597","SLE_GSE263931_GSM8207599","SLE_GSE263931_GSM8207601","SLE_GSE263931_GSM8207603","SLE_GSE263931_GSM8207605","SLE_GSE263931_GSM8207607","SLE_GSE224198_GSM7017326","SLE_GSE224198_GSM7017329","SLE_GSE224198_GSM7017331","SLE_GSE224198_GSM7017334","JIA_GSE205095_GSM6205132","JIA_GSE205095_GSM6205136","JIA_GSE205095_GSM6205138","JIA_GSE205095_GSM6205140","JIA_GSE205095_GSM6205142","JIA_GSE205095_GSM6205144","KD_GSE168732_GSM5160417","KD_GSE168732_GSM5160420","KD_GSE168732_GSM5160422","KD_GSE168732_GSM5160424","KD_GSE168732_GSM5160427","KD_GSE168732_GSM5160430","STING_GSE226598_GSM7079986","STING_GSE226598_GSM7079989","STING_GSE226598_GSM7079991","PLCG2_Figshare","OTULIN_GSE199445_GSM5973147","OTULIN_GSE199445_GSM5973148","OTULIN_GSE199445_GSM5973149")
metainfo<-fread("metainfo.txt")
fold2<-getwd()


Quality_check_debug <- function(samples, fold2, outdir = NULL, pdf_file = NULL) {
  
  message("=====QC check start====")
  start_time <- Sys.time()  
  
  if (!dir.exists(outdir)) dir.create(outdir)
  res_list <- lapply(samples, function(x) {
    message("\n--- Processing sample: ", x, " ---")
    
    s <- CreateSeuratObject(counts = Read10X(file.path(fold2, x)))
    s[["percent.mt"]] <- PercentageFeatureSet(s, pattern = "^MT-")
    s[["percent.ribo"]] <- PercentageFeatureSet(s, pattern = "^RPS|^RPL")
    s[["percent.platelet"]] <- PercentageFeatureSet(s, features = c("PF4","PPBP","ITGA2B","GNG11","TUBB1"))
    s[["percent.rbc"]] <- PercentageFeatureSet(s, features = c("HBB", "HBD", "HBA1", "HBA2"))
   
    meta <- s@meta.data
    
    print(head(meta))
    
    # check if there exist NA in metadata
    na_counts <- colSums(is.na(meta))
    if (any(na_counts > 0)) {
      warning("Sample ", x, " has NA in the following columns:")
      print(na_counts[na_counts > 0])
    }
    
    # check if there exit one column only containing one constant
    constant_cols <- sapply(meta, function(c) length(unique(c)) == 1)
    if (any(constant_cols)) {
      message("Sample ", x, " has constant value in columns: ", paste(names(meta)[constant_cols], collapse=", "))
    }
    
    # Vioplots
    p <- VlnPlot(s,
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo","percent.platelet", "percent.rbc"),
                 ncol = 4,
                 pt.size = 0) & theme(axis.title.x = element_blank(),
                                      axis.text.x = element_blank(),
                                      axis.ticks.x = element_blank())
    
    p<-p+ plot_annotation(
      title = paste(x,as.character(metainfo[match(x,metainfo$dataset),'chemistry']),sep=", "),  
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
      ))
    
    
    # combine to a summary DF
    df <- data.frame(
      sample = x,
      n_cells = ncol(s),
      median_nFeature = median(meta$nFeature_RNA, na.rm = TRUE),
      median_nCount = median(meta$nCount_RNA, na.rm = TRUE),
      median_mt = median(meta$percent.mt, na.rm = TRUE),
      median_ribo = median(meta$percent.ribo, na.rm = TRUE),
      median_platelet = median(meta$percent.platelet, na.rm = TRUE),
      median_rbc = median(meta$percent.rbc, na.rm = TRUE),
      mean_nFeature = mean(meta$nFeature_RNA, na.rm = TRUE),
      mean_nCount = mean(meta$nCount_RNA, na.rm = TRUE),
      mean_mt = mean(meta$percent.mt, na.rm = TRUE),
      mean_ribo = mean(meta$percent.ribo, na.rm = TRUE),
      mean_platelet = mean(meta$percent.platelet, na.rm = TRUE),
      mean_rbc = mean(meta$percent.rbc, na.rm = TRUE),
      min_nFeature = min(meta$nFeature_RNA, na.rm = TRUE),
      min_nCount = min(meta$nCount_RNA, na.rm = TRUE),
      min_mt = min(meta$percent.mt, na.rm = TRUE),
      min_ribo = min(meta$percent.ribo, na.rm = TRUE),
      min_platelet = min(meta$percent.platelet, na.rm = TRUE),
      min_rbc = min(meta$percent.rbc, na.rm = TRUE),
      max_nFeature = max(meta$nFeature_RNA, na.rm = TRUE),
      max_nCount = max(meta$nCount_RNA, na.rm = TRUE),
      max_mt = max(meta$percent.mt, na.rm = TRUE),
      max_ribo = max(meta$percent.ribo, na.rm = TRUE),
      max_platelet = max(meta$percent.platelet, na.rm = TRUE),
      max_rbc = max(meta$percent.rbc, na.rm = TRUE),
      p01_nFeature = quantile(meta$nFeature_RNA, 0.01, na.rm = TRUE),
      p99_nFeature = quantile(meta$nFeature_RNA, 0.99, na.rm = TRUE),
      p01_mt = quantile(meta$percent.mt, 0.01, na.rm = TRUE),
      p99_mt = quantile(meta$percent.mt, 0.99, na.rm = TRUE),
      p01_ribo = quantile(meta$percent.ribo, 0.01, na.rm = TRUE),
      p99_ribo = quantile(meta$percent.ribo, 0.99, na.rm = TRUE),
      p01_platelet = quantile(meta$percent.platelet, 0.01, na.rm = TRUE),
      p99_platelet = quantile(meta$percent.platelet, 0.99, na.rm = TRUE),
      p01_rbc = quantile(meta$percent.rbc, 0.01, na.rm = TRUE),
      p99_rbc = quantile(meta$percent.rbc, 0.99, na.rm = TRUE)
    )
    
    
    return(list(plot = p, summary = df))
  })
  
  # combine summary
  summary_df <- do.call(rbind, lapply(res_list, function(x) x$summary))
  summary_df$chemistry<-as.character(as.matrix(metainfo[match(summary_df$sample,metainfo$dataset),'chemistry']))
  summary_df<- as.data.frame(lapply(summary_df, as.character), stringsAsFactors = FALSE)
  if (!is.null(outdir)) {
    write.xlsx(summary_df, file = file.path(outdir, "scRNA_rawQCdata_debug.xlsx"))
  }
  
  # combine raw qc violin plots for all samples
  plot_list<-lapply(res_list,function(x){
    ggdraw()+draw_plot(x$plot)})
  
  pdf(file.path(outdir, pdf_file), width = 18, height = 20)
  n_per_page <- 8  # ncol 2 nrow 4
  for (i in seq(1, length(plot_list), by = n_per_page)) {
    
    page_plot <- plot_grid(
      plotlist = plot_list[i:min(i + n_per_page - 1, length(plot_list))],
      ncol = 2
    )
    
    print(page_plot)
  }
  dev.off()
  
  message("=====QC check done==== ", normalizePath(outdir))
  end_time <- Sys.time()
  message("Total runtime: ", round(difftime(end_time, start_time, units = "mins"), 2), " minutes")
  
  return(summary_df)
}
summary_debug <- Quality_check_debug(samples, fold2, outdir = fold2, pdf_file = "quality_summary_debug.pdf")