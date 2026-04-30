library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(ggplot2)


samples=c("137ACP5","84ADA2","159ADA2","163COPA","88GNAI2","78IL1R1","79IL1R1_2","156ISG15","157LACC1","90OGFRL1","36PLCG1","37PLCG2","108PLD4","109PLD4","110PLD4","155PSMD12","138RIPK1","139RIPK1","140RIPK1","141RIPK1","152RNF31","149SOCS1","162STING","69TBXAS1","70TBXAS1","89TLR7","82TNFAIP3","83TNFAIP3","154TRAF3","151USP18","134ELF4","115KRAS","81NLRC4","91TLR8","148UBA1","1IRAK2","3IRAK2","165TREX1","168TREX1","169TREX1","175IFIH1","176IFIH1","177IFIH1","178ISG15","179STAT4","180SLC7A7","181NOD2","182LPIN2","12HC","17HC","22HC","75HC","76HC","77HC","158HC","160HC","161HC","9HC","HC_GSE139324_GSM4138162","HC_GSE139324_GSM4138163","HC_GSE139324_GSM4138164","HC_GSE139324_GSM4138165","HC_GSE139324_GSM4138166","HC_GSE139324_GSM4138167","HC_10XGenomics10K","HC_GSE156989_GSM4749756","HC_GSE156989_GSM4749762","HC_GSE156989_GSM4749768","HC_GSE199445_GSM5973143","HC_GSE199445_GSM5973144","HC_GSE199445_GSM5973145","HC_GSE199445_GSM5973146","HC_GSE168732_GSM5160432","HC_GSE168732_GSM5160434","HC_GSE168732_GSM5160435","SLE_GSE156989_GSM4749774","SLE_GSE156989_GSM4749779","SLE_GSE156989_GSM4749784","SLE_GSE142016_GSM4217718","SLE_GSE142016_GSM4217719","SLE_GSE142016_GSM4217720","SLE_GSE263931_GSM8207595","SLE_GSE263931_GSM8207597","SLE_GSE263931_GSM8207599","SLE_GSE263931_GSM8207601","SLE_GSE263931_GSM8207605","SLE_GSE263931_GSM8207607","SLE_GSE224198_GSM7017326","SLE_GSE224198_GSM7017329","SLE_GSE224198_GSM7017331","SLE_GSE224198_GSM7017334","JIA_GSE205095_GSM6205132","JIA_GSE205095_GSM6205136","JIA_GSE205095_GSM6205138","JIA_GSE205095_GSM6205140","JIA_GSE205095_GSM6205142","JIA_GSE205095_GSM6205144","KD_GSE168732_GSM5160417","KD_GSE168732_GSM5160420","KD_GSE168732_GSM5160422","KD_GSE168732_GSM5160424","KD_GSE168732_GSM5160427","KD_GSE168732_GSM5160430","STING_GSE226598_GSM7079986","STING_GSE226598_GSM7079989","STING_GSE226598_GSM7079991","PLCG2_Figshare","OTULIN_GSE199445_GSM5973147","OTULIN_GSE199445_GSM5973148","OTULIN_GSE199445_GSM5973149")
input="/public/home/yuxiaomingroup/msy/work/AID/juzhen/clean/result/"

#load h5Seurat files to combine a large Seurat object
process_HVG_Lisi<-function(input,samples){
  #sce_list<-lapply(samples,function(x){file=paste0(input,x,"Redl.h5Seurat")
  #s<-LoadH5Seurat(file)})
  
  #saveRDS(sce_list,file="reDl_scellist.RDS")
  
  sce_list<-readRDS("reDl_scellist.RDS")
  # merge cells and genes by taking the union
  #pbmc <- merge(x = sce_list[[1]],y = sce_list[-1],add.cell.ids =samples,project = "PBMC")
  
  #pbmc <- DietSeurat(pbmc,counts = TRUE,data = TRUE,scale.data = FALSE)
  #SaveH5Seurat(pbmc, filename = "pbmcmerged.h5Seurat")
  pbmc<-LoadH5Seurat("pbmcmerged.h5Seurat")
  
  ## (1) HVG genes was selected by top 3000 rank across the union hgvs derived from each dataset
  pbmc1 <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  hvg_list <- lapply(sce_list, VariableFeatures)# select shared HVG across list elements
  hvgs <- Reduce("union",hvg_list)
  hvgs <- hvgs[hvgs %in% rownames(pbmc1)]
  length(hvgs)
  gene_count <- sapply(hvgs, function(g){sum(sapply(hvg_list, function(hvgs_sample) g %in% hvgs_sample))})
  threshold <- sort(gene_count, decreasing = TRUE)[3000]
  hvgs_top3000 <- names(gene_count[gene_count >= threshold])
  gene_count_sorted <- sort(gene_count, decreasing = TRUE)
  df <- data.frame(rank = seq_along(gene_count_sorted),count = gene_count_sorted)
  
  hvgCount_FreqBar<-ggplot(df, aes(x = rank, y = count)) +
    geom_line(color = "#782123", linewidth = 0.5) +
    geom_vline(xintercept = 3000, linetype = "dashed",linewidth=0.3) +
    labs(x = "Gene rank",y = "Number of datasets",title = "")+theme(
      axis.text = element_text(color = "black"))+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)),
                       limits = c(0, 110))+theme_classic()+theme(
                         axis.text = element_text(colour = "black"),
                         axis.ticks = element_line(colour = "black"))+  annotate(
                           "text",x = 7400,y = 110,
                           label = paste0("top 3000 rank: \n", length(hvgs_top3000), " hvgs"),
                           color = "black",size = 3) 
  
  ##(2) HVG genes was derived by find variable genes for merged seurat, run one time
  merge_hvgs<-VariableFeatures(FindVariableFeatures(pbmc1, selection.method = "vst", nfeatures = 3000))
  
  ggsave(
    "hvgCount_FreqBar.pdf",
    plot = hvgCount_FreqBar,
    width = 5,
    height = 4.5,
    units = "in"
  )
  saveRDS(merge_hvgs, file = "merge_hvgs.rds")
  saveRDS(hvgs_top3000, file = "hvgs_top3000.rds")
  saveRDS(gene_count, file = "gene_count.rds")
  saveRDS(hvg_list, file = "hvg_list.rds")
}
process_HVG_Lisi(input,samples)