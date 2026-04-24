
library(Seurat)
library(SeuratDisk)
library(openxlsx)

samples=c("137ACP5","84ADA2","159ADA2","163COPA","88GNAI2","78IL1R1","79IL1R1_2","156ISG15","157LACC1","90OGFRL1","36PLCG1","37PLCG2","108PLD4","109PLD4","110PLD4","155PSMD12","138RIPK1","139RIPK1","140RIPK1","141RIPK1","152RNF31","149SOCS1","162STING","69TBXAS1","70TBXAS1","89TLR7","82TNFAIP3","83TNFAIP3","154TRAF3","151USP18","134ELF4","115KRAS","81NLRC4","91TLR8","148UBA1","1IRAK2","3IRAK2","165TREX1","168TREX1","169TREX1","175IFIH1","176IFIH1","177IFIH1","178ISG15","179STAT4","180SLC7A7","181NOD2","12HC","17HC","22HC","75HC","76HC","77HC","158HC","160HC","161HC","9HC","HC_GSE139324_GSM4138162","HC_GSE139324_GSM4138163","HC_GSE139324_GSM4138164","HC_GSE139324_GSM4138165","HC_GSE139324_GSM4138166","HC_GSE139324_GSM4138167","HC_10XGenomics10K","HC_GSE156989_GSM4749756","HC_GSE156989_GSM4749762","HC_GSE156989_GSM4749768","HC_GSE199445_GSM5973143","HC_GSE199445_GSM5973144","HC_GSE199445_GSM5973145","HC_GSE199445_GSM5973146","HC_GSE168732_GSM5160432","HC_GSE168732_GSM5160434","HC_GSE168732_GSM5160435","SLE_GSE156989_GSM4749774","SLE_GSE156989_GSM4749779","SLE_GSE156989_GSM4749784","SLE_GSE142016_GSM4217718","SLE_GSE142016_GSM4217719","SLE_GSE142016_GSM4217720","SLE_GSE263931_GSM8207595","SLE_GSE263931_GSM8207597","SLE_GSE263931_GSM8207599","SLE_GSE263931_GSM8207601","SLE_GSE263931_GSM8207603","SLE_GSE263931_GSM8207605","SLE_GSE263931_GSM8207607","SLE_GSE224198_GSM7017326","SLE_GSE224198_GSM7017329","SLE_GSE224198_GSM7017331","SLE_GSE224198_GSM7017334","JIA_GSE205095_GSM6205132","JIA_GSE205095_GSM6205136","JIA_GSE205095_GSM6205138","JIA_GSE205095_GSM6205140","JIA_GSE205095_GSM6205142","JIA_GSE205095_GSM6205144","KD_GSE168732_GSM5160417","KD_GSE168732_GSM5160420","KD_GSE168732_GSM5160422","KD_GSE168732_GSM5160424","KD_GSE168732_GSM5160427","KD_GSE168732_GSM5160430","STING_GSE226598_GSM7079986","STING_GSE226598_GSM7079989","STING_GSE226598_GSM7079991","PLCG2_Figshare","OTULIN_GSE199445_GSM5973147","OTULIN_GSE199445_GSM5973148","OTULIN_GSE199445_GSM5973149")
input="E:\\AID cohort\\code\\result\\"

filter2Seurat<-function(samples,input){
  
message("=====saving transform start====")
start_time <- Sys.time() 

douletDF<-lapply(samples,function(x){

s <- CreateSeuratObject(Read10X(paste0(input,"filter",x)))
obs<-read.xlsx(paste0(input,"filter",x,"/cell_metadata.xlsx"))
var<-read.xlsx(paste0(input,"filter",x,"/gene_metadata.xlsx"))

s@meta.data$nCount_RNA<-obs[match(rownames(s@meta.data),obs$X1),'n_counts']
s@meta.data$nFeature_RNA<-obs[match(rownames(s@meta.data),obs$X1),'n_genes']
s@meta.data<-cbind(s@meta.data,
                   obs[match(rownames(s@meta.data),obs$X1),
                       c("samples","percent_ribo","percent_mt","log10GenesPerUMI","doublet_score","predicted_doublet")])
s@assays$RNA@var.features<-as.character(
  var[var$highly_variable==TRUE,'X1'])

s1 <- CreateSeuratObject(Read10X(paste0(input,"Counts",x)))
s@assays$RNA@counts<-(s1@assays$RNA@counts)
rm(s1);gc()
scores<-s@meta.data[,c('samples',"doublet_score")]

setwd(input)
SaveH5Seurat(DietSeurat(s),filename = paste0(x,"Redl.h5Seurat"))
message(paste0("====",x,"has been save to h5Seurat"))
rm(s);gc()
return(scores)
})
message("=====QC done====") # 1 as TRUE as doublet
end_time <- Sys.time()  
message("Total runtime: ", round(difftime(end_time, start_time, units = "mins"), 2), " minutes")
return(douletDF)
}

results<-filter2Seurat(samples,input)

