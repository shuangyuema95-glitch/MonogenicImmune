library(Seurat)
library(SeuratDisk)
pbmc1<-readRDS("/public/home/yuxiaomingroup/msy/work/AID/juzhen/clean/result/50/PBMC1.rds")
pbmc1$Level1<-as.character(pbmc1$Level1)
DefaultAssay(pbmc1) <- "RNA"
pbmc1[["RNA"]]@scale.data <- matrix()
sum(is.na(GetAssayData(pbmc1, slot = "data")))
SaveH5Seurat(pbmc1, filename = "716_pbmc1.h5Seurat", overwrite = TRUE)
Convert("716_pbmc1.h5Seurat", dest = "h5ad", overwrite = TRUE)