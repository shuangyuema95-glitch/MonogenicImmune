#' Sublineage Secondary Clustering Pipeline
#' @title run_single_lineage: Recluster Single Major Cell Subset with Multiple Resolutions + Harmony UMAP Visualization
#' @description Extract specified Level1 cell subset from full integrated Seurat object, re-run PCA/Harmony/UMAP/clustering under 3 resolutions (0.05,0.08,0.1). Generate cluster UMAP & sample batch UMAP plots, save subsetted Seurat rds for each resolution. Return cell count table of level2 clusters at resolution 0.1.
#' Step1: Load full seurat object, rename original harmony clusters to Level1 cell type annotation
#' Step2: Subset target lineage, clear original reduction slots to avoid interference
#' Step3: QC filter: skip clustering if total cell number < 50
#' Step4: Standard single-cell workflow: HVG selection, regress out cell cycle & mt genes, PCA + Harmony batch correction
#' Step5: KNN graph construction, multi-resolution clustering & harmony-based UMAP embedding
#' Step6: Plot cluster UMAP with centroid labels, sample batch UMAP without legend, export 4x4 pdf figures
#' Step7: Save subsetted Seurat object annotated with level2 cluster metadata for each resolution
#' Step8: Load rds of res=0.1, output cluster cell count table as final return value
#' @param full_cell_name Character, target Level1 cell type to subset from full object
#' @param short_tag Character, prefix tag for all output figure & rds filenames
#' @returnType Table
#' @return Integer table showing cell counts per level2 cluster under resolution 0.1
#' @author MSY
run_single_lineage <- function(full_cell_name,short_tag){
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  res_list <- c(0.05,0.08,0.1)
  for(this_res in res_list){
    message(paste0("==== Start running resolution = ", this_res, " ===="))
    pbmc1 <- readRDS(seu_path)
    Idents(pbmc1) <- "harmony_clusters"
    pbmc1 <- RenameIdents(pbmc1,'0'='Naive CD4 T','1'='Non-naive CD8 T','2'='B cell','3'='Non-naive CD4 T','4'='Naive CD8 T','5'='Monocyte','6'='NK','7'='Monocyte','8'='LDG','9'='Monocyte','10'='Platelet','11'='Monocyte','12'='Naive CD4 T','13'='Erythrocyte','14'='Plasma','15'='Monocyte','16'='DC','17'='UTC','18'='Monocyte','19'='Basophil','20'='pDC','21'='Monocyte','22'='B cell','23'='Monocyte','24'='B cell','25'='Monocyte','26'='Erythrocyte','27'='Non-naive CD4 T','28'='DC','29'='DC','30'='Basophil','31'='Monocyte','32'='Monocyte','33'='Monocyte')
    pbmc1$Level1 <- Idents(pbmc1)
    sub_obj <- subset(pbmc1,Level1==full_cell_name)
    sub_obj@reductions <- list()
    rm(pbmc1);gc()
    if(ncol(sub_obj)<50){message(sprintf("[%s] Cell count <50, skip clustering",full_cell_name));return(NULL)}
    pcaNumber <- 30
    sub_obj <- FindVariableFeatures(sub_obj,selection.method="vst",nfeatures=2000)
    sub_obj <- ScaleData(sub_obj,features=VariableFeatures(sub_obj),vars.to.regress=c("percent_mt","S.Score","G2M.Score","nCount_RNA"))
    sub_obj <- RunPCA(sub_obj,features=VariableFeatures(sub_obj),npcs=pcaNumber)
    sub_obj <- RunHarmony(sub_obj,group.by.vars="samples",dims.use=1:pcaNumber,theta=2)
    sub_obj <- FindNeighbors(sub_obj,dims=1:pcaNumber,reduction="harmony",k.param=30)
    sub_obj <- FindClusters(sub_obj,resolution=this_res)
    sub_obj <- RunUMAP(sub_obj,reduction="harmony",dims=1:pcaNumber,reduction.name="harmony_umap",key="harUMAP_")
    plot_df <- as.data.frame(Embeddings(sub_obj,"harmony_umap"))
    colnames(plot_df) <- c("UMAP_1","UMAP_2")
    plot_df$Cluster <- as.factor(sub_obj$seurat_clusters)
    plot_df$samples <- sub_obj$samples
    centroid_df <- plot_df %>% group_by(Cluster) %>% summarise(UMAP_1=mean(UMAP_1),UMAP_2=mean(UMAP_2),.groups="drop")
    p_cluster <- ggplot(plot_df,aes(UMAP_1,UMAP_2,color=Cluster))+geom_point(size=0.15,shape=16,stroke=0)+scale_color_manual(values=sub_cluster_colors)+geom_text_repel(data=centroid_df,aes(label=Cluster),size=3,fontface="bold",color="black",max.overlaps=Inf,show.legend=FALSE)+theme_classic()+theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.position="none",aspect.ratio=1)+ggtitle(paste0(full_cell_name," | Cluster, res=", this_res, " N = ",ncol(sub_obj)))
    ggsave(file.path(out_dir,paste0(short_tag,"_res",this_res,"_level2_cluster_umap.pdf")),p_cluster,width=4,height=4)
    p_sample <- ggplot(plot_df,aes(UMAP_1,UMAP_2,color=samples))+geom_point(size=0.15,shape=16,stroke=0)+theme_classic()+theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.position="none",aspect.ratio=1)+guides(color="none")+ggtitle(paste0(full_cell_name," | Sample Batch Check, res=", this_res))
    ggsave(file.path(out_dir,paste0(short_tag,"_res",this_res,"_level2_sample_umap.pdf")),p_sample,width=4,height=4)
    sub_obj$level2_cluster <- sub_obj$seurat_clusters
    saveRDS(sub_obj,file.path(out_dir,paste0(short_tag,"_res",this_res,"_level2_seu.rds")))
    message(paste0("==== Finish resolution = ", this_res," ====\n"))
  }
  final_res0.1 <- readRDS(file.path(out_dir,paste0(short_tag,"_res0.1_level2_seu.rds")))
  return(table(final_res0.1$level2_cluster))
}

library(Seurat)
library(tidyverse)
library(ggrepel)
library(cowplot)
library(harmony)

cellcase <- "col1"
fexcase <- "col2"
out_dir <- "./sublineage_result"
seu_path <- "/public/home/yuxiaomingroup/msy/work/AID/juzhen/clean/result/50/50_t3k_full_backup.rds"
sub_cluster_colors <- c("#847AB3","#F6E1EE","#B3928B","#C8BBD5","#C95968","#F49D5C","#B3D1E7","#96AF95","#003D81","#C6D199","#74AB9B","#D6D78D","#0B71AB","#8C2522","#FFEFC1")

run_single_lineage(full_cell_name=cellcase,short_tag=fexcase)