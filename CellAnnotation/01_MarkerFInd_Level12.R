#' Marker Identification Pipeline for Full Cohort Seurat Object
#' @title RunAllMarkers_Cohort: FindAllMarkers + FindConservedMarkers + COSG
#' @description Perform three marker detection modules on integrated seurat object, parameters fixed as original setting.
#' Step1: FindAllMarkers with only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25; export rds + xlsx.
#' Step2: FindConservedMarkers grouped by sample batch, add gene annotation and cluster_id column, export rds + xlsx.
#' Step3: COSG cluster-specific marker calculation, fixed parameters, export raw result rds and top gene excel.
#' All intermediate results saved separately to avoid full loss, return combined marker list object.
#' @param seurat_obj Integrated Seurat object
#' @param outdir Output directory for marker tables
#' @param annotations Gene annotation dataframe with gene_name, description columns
#' @param cluster_col Metadata column storing cluster id
#' @param batch_var Metadata column for sample grouping in conserved marker test
#' @returnType List
#' @return List containing all_markers dataframe, conserved_marker full dataframe, conserved marker gene list, raw COSG output
#' @author MSY


library(Seurat)
library(SeuratObject)
library(ggplot2)
library(tidyverse)
library(purrr)
library(COSG)
library(openxlsx)
ROOT <- "/public/home/yuxiaomingroup/msy/work/AID/juzhen/clean/result/"
pcaNumber <- 50
sufix <- "t3"
targetfold <- file.path(ROOT, as.character(pcaNumber))
outdir <- file.path(targetfold, "markers")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
load(file.path(ROOT, "Marker_Annotation.Rdata"))

options(future.globals.maxSize = 100 * 1024^3)
set.seed(123)

RunAllMarkers_Cohort <- function(seurat_obj, outdir, annotations, cluster_col = "harmony_clusters", batch_var = "samples") {
  pbmc1 <- seurat_obj
  Idents(pbmc1) <- cluster_col
  
  message("==== Step1 FindAllMarkers ====")
  markers <- FindAllMarkers(pbmc1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  saveRDS(markers, file.path(outdir, paste0(pcaNumber, "_", sufix, "_AllMarkers.rds")))
  write.xlsx(markers, file.path(outdir, paste0(pcaNumber, "_", sufix, "_AllMarkers.xlsx")), rowNames = FALSE)
  
  message("==== Step2 FindConservedMarkers across samples ====")
  get_conserved <- function(cluster) {
    cm <- FindConservedMarkers(pbmc1, ident.1 = cluster, grouping.var = batch_var, ident.2 = NULL, only.pos = TRUE) %>%
      rownames_to_column(var = "gene") %>%
      left_join(unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name")) %>%
      mutate(cluster_id = as.character(cluster))
    return(cm)
  }
  cluster_vec <- sort(unique(pbmc1[[cluster_col, drop = TRUE]]))
  con_markers_df <- map_dfr(cluster_vec, get_conserved)
  con_markers_list <- split(con_markers_df$gene, f = con_markers_df$cluster_id)
  saveRDS(con_markers_list, file.path(outdir, paste0(pcaNumber, "_", sufix, "_ConservedMarkerList.rds")))
  write.xlsx(con_markers_df, file.path(outdir, paste0(pcaNumber, "_", sufix, "_ConservedMarkerFull.xlsx")), rowNames = FALSE)
  
  message("==== Step3 COSG cluster marker ====")
  cosg_res <- cosg(pbmc1, groups = "all", assay = "RNA", slot = "data", mu = 10, n_genes_user = 400, remove_lowly_expressed = TRUE, expressed_pct = 0.1)
  cosg_df <- as.data.frame(cosg_res$names)
  saveRDS(cosg_res, file.path(outdir, paste0(pcaNumber, "_", sufix, "_COSG_Result.rds")))
  write.xlsx(cosg_df, file.path(outdir, paste0(pcaNumber, "_", sufix, "_COSG_TopGenes.xlsx")), rowNames = FALSE)
  
  return(list(all_markers = markers, conserved_marker_df = con_markers_df, conserved_marker_list = con_markers_list, cosg_result = cosg_res))
}

seu_file <- file.path(targetfold, paste0(pcaNumber, "_", sufix, "_allReduction_full.h5Seurat"))
obj <- LoadH5Seurat(seu_file)

marker_result <- RunAllMarkers_Cohort(seurat_obj = obj, outdir = outdir, annotations = annotations, cluster_col = "harmony_clusters", batch_var = "samples")

cat("Marker analysis finished!\nOutput dir:", outdir, "\n")
cat("Total clusters:", length(unique(obj[[cluster_col, drop=T]])), "\n")