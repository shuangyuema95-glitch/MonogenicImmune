#################
################# Testing correlation with CellTypist
library(Seurat)
library(pheatmap)
library(tidyverse)
library(viridis)

# ====================== 1. Define all marker gene sets ======================
markerlist<-list(
  "Naive CD4 T" = c('LEF1', 'CCR7', 'TCF7',"MAL","IL7R","CD4"),
  "Non-naive CD8 T"= c("CD8B","GZMH","GZMA","CCL5","FCRL6"),
  "B cell"=c("CD19","MS4A1", "CD79A", "CD79B", "BANK1", "CD22"),
  "Non-naive CD4 T"=c("IL32","AQP3","GATA3","RORA","ICOS"),
  "Naive CD8 T"=c("CD8B","CD8A","CCR7","LEF1","TCF7","LDLRAP1"),
  "Monocyte"=c("CD14", "S100A8", "S100A9","FCGR3A","LST1","LILRB2"),
  "NK"=c("NKG7","KLRF1","KIR2DL3","KIR3DL1","NCAM1"),
  "LDG"=c("FCGR3B","CSF3R","SRGN","IFITM2","CMTM2","FUT4"),
  "Platelet"=c("PPBP","PF4","TUBB1"),
  "Erythrocyte"=c("HBA1","HBA2","HBB"),
  "Plasma" =c ("JCHAIN","MZB1","IGHA1","XBP1"),
  "DC"=c("FCER1A","CD1C","HLA-DRA","HLA-DQA1"),
  "UTC"=c("GZMK","TRGC2","TRAC","CCL5","TIGIT","CXCR6"),
  "Basophil" = c("CCR3","MS4A3","CLC","ENPP3"),
  "pDC"=c("PLD4","IL3RA","LILRA4","TCF4")
)

B_markerlist<-list(
  B01_Naive_BACH2=c("BACH2","FCRL1","IGHD","SNX29","IGHM","MS4A1"),
  B02_Transitional_TCL1A=c("TCL1A","VPREB3","IGLL5","FCRLA","HVCN1","FCER2"),
  B03_MemorySwitch_SNED1=c("SNED1","TNFRSF13B", "ARHGAP24","COCH","AIM2","CRIP2","BLK"),
  "B04_CD27+BCR+"=c("CD27","BCR","XBP1","IGLC1","CD52","IGHA1","IGHA2"),
  B05_IFNresponse_MX1=c("MX1","IFI44L","IFIT1","OAS1","SAMD9L"),
  B06_MemorySwitch_IGHA1=c("IGHA1","IGHA2","IGKC")
)

DC_markerlist<-list(
  DC01_cDC2_CD1C=c("CD1C","FCER1A","CLEC10A","CD1E"),
  DC02_DC_TRIM33=c("TRIM33","MAML3","PLXDC2"),
  DC03_aDC_CCR7=c("CCR7","CD40","CCL22","CD83","LAMP3"),
  DC04_DC_CCL5=c("CCL5","CTSD","CLEC2D"),
  DC05_cDC2_LY96=c("LY96","CSTA","FAM26F"),
  DC06_cDC1_XCR1=c("CADM1","XCR1","CLEC9A","IRF8","BATF3"),
  DC07_mDC_ISG15=c("ISG15","LILRB1","FCGR2B","CD14","CD68"),
  DC08_DC_C1QTNF4 =c("C1QTNF4","LAIR1","IGFBP7"),
  DC09_mDC_VCAN=c("VCAN","LYZ","CD1C","LGALS3","CD14"))

NK_markerlist<-list(
  NK01_CD56dimCD16_KLRC2=c("NCAM1","FCGR3A","PITPNC1","KLRC2"),
  NK02_CD16high_PRF1=c("FCGR3A","NKG7","PRF1"),
  NK03_CD56high_GZMK=c("NCAM1","GZMK","XCL1"),
  NK04_IFNresponse_MX1=c("MX1","OAS1","IFIT1","CXCL10","KLRF1"))

Plasma_markerlist<-list(
  Plasma01_JCHAIN=c("JCHAIN","IGHG1","MZB1","IGLC3"),
  Plasma02_CXCR4=c("CXCR4","IGHA1","IGLC2"),
  Plasma03_Cycling_MKI67=c("TOP2A","MKI67"),
  Plasma04_XBP1=c("XBP1","IGHG4","IGHG2")
)

NnCD4T_markerlist<-list(
  NnCD4T01_Memory_RORA=c("RORA","CREM","ZNF331","BACH2","FOXP1"),
  NnCD4T02_CM_FOS=c("FOS","CD27","TCF7","COTL1","CRIP1"),
  NnCD4T03_Treg_IKZF2=c("IKZF2","IL2RA","FOXP3"),
  NnCD4T04_CM_IFN_HERC5=c("HERC5","IFIH1","DDX58"),
  NnCD4T05_EMRA_GNLY=c("GNLY","GZMA","GZMK","CCL5","NKG7"),
  NnCD4T06_CM_IFN_CISH=c("CD27","CCR7","CISH","RTP4","IFIT1"))

NnCD8T_markerlist<-list(
  NnCD8T01_Effector_ZEB2=c("ZEB2","EOMES","IFNG","TBX21"),
  NnCD8T02_Cytotoxic_PRF1=c("PRF1","FGFBP2","GZMH","GZMB","GNLY","NKG7"),
  NnCD8T03_CM_NELL2=c("NELL2","RCAN3","FOXP1","ZNF331"),
  NnCD8T04_EM_CXCR6=c("RORA","KLRB1","CXCR6","GZMK","IL7R"),
  NnCD8T05_EM_CMC1=c("CMC1","CTSW","CCL5"),
  NnCD8T06_CM_IFN_MX1=c("MX1","OAS1","HERC6","SELL","CD27")
)

UTC_markerlist<-list(
  UTC01_MAIT_GZMK=c("GZMK","CXCR3","RORA"),
  UTC02_gdTV2Vd9=c("BTN3A1"),
  UTC03_gdTV1=c("TRDC","TRGC1","KLRB1"))

Mono_markerlist<-list(
  Mono01_Classical_FOSB=c("FOSB","CEBPD","CD14"),
  Mono02_Inflammatory=c("CXCL8","IL1B","NLRP3","ICAM1","CD44"),
  Mono03_Nonclassical=c("FCGR3A","MS4A7","LST1","CSF3R","CD68"),
  Mono04_Classical_RETN=c("MYDGF","RETN","FCHSD2","CD14"),
  Mono05_Classical_S100A9=c("S100A9","S100A8","S100A12"),
  Mono06_Classical_TLN1=c("TLN1","MYL9","PGRMC1"),
  Mono07_IFNresponse=c("S100A8","LYZ","IFIT1","RSAD2","LY6E","IFI44L"))

# Combine all markers and deduplicate
all_marker_genes <- unique(c(
  unlist(markerlist), unlist(B_markerlist), unlist(DC_markerlist),
  unlist(NK_markerlist), unlist(Plasma_markerlist), unlist(NnCD4T_markerlist),
  unlist(NnCD8T_markerlist), unlist(UTC_markerlist), unlist(Mono_markerlist)
))

# Function to filter existing genes
gene_filter <- function(mat, gene_vec){
  inter_gene <- intersect(gene_vec, rownames(mat))
  mat[inter_gene, , drop = FALSE]
}

# ====================== 2. Load Seurat object and merge celltypist label ======================
pbmc1<-readRDS("/public/home/yuxiaomingroup/msy/work/AID/juzhen/clean/result/50/PBMC1.rds")
meta_cp<-read.csv("/public/home/yuxiaomingroup/msy/work/AID/juzhen/clean/result/50/pbmc1_celltypist_metadata.csv")

pbmc1$celltypist_label <- meta_cp$celltypist_label

# ====================== 3. Calculate full average expression matrices ======================
avg_lv2_full <- AverageExpression(pbmc1, group.by = "Level2", assay = "RNA")$RNA
avg_lv2_full <- as.matrix(avg_lv2_full)

avg_ct_full <- AverageExpression(pbmc1, group.by = "celltypist_label", assay = "RNA")$RNA
avg_ct_full <- as.matrix(avg_ct_full)

# ====================== 4. Save raw average matrices and free huge Seurat memory ======================
saveRDS(avg_lv2_full, "avg_lv2_full.rds")
saveRDS(avg_ct_full, "avg_ct_full.rds")
rm(pbmc1, meta_cp)
gc()

# ====================== 5. Load matrices and filter cell subtypes ======================
avg_lv2_full <- readRDS("avg_lv2_full.rds")
avg_ct_full <- readRDS("avg_ct_full.rds")

# Filter celltypist blacklist
black_ct <- c("Mast cells","Intermediate macrophages","HSC/MPP","Alveolar macrophages")
ct_keep <- setdiff(colnames(avg_ct_full), black_ct)
avg_ct_filter <- avg_ct_full[, ct_keep, drop = F]
rm(avg_ct_full)
gc()

# Filter Level2 blacklist: Basophil, LDG
lv2_black <- c("Basophil", "LDG")
lv2_keep <- setdiff(colnames(avg_lv2_full), lv2_black)
avg_lv2_filter <- avg_lv2_full[, lv2_keep, drop = F]
rm(avg_lv2_full)
gc()

# ====================== 6. Subset matrices to marker genes only ======================
avg_lv2_marker <- gene_filter(avg_lv2_filter, all_marker_genes)
avg_ct_marker  <- gene_filter(avg_ct_filter, all_marker_genes)
rm(avg_lv2_filter, avg_ct_filter)
gc()

# ====================== 7. Compute correlation matrix ======================
expr_bind <- cbind(avg_lv2_marker, avg_ct_marker)
cor_all <- cor(expr_bind)

# Extract plotting matrix: row = Level2, column = celltypist label
cor_plot <- cor_all[colnames(avg_lv2_marker), colnames(avg_ct_marker)]

# Save cor_plot for local plotting adjustment later
saveRDS(cor_plot, "cor_plot_matrix.rds")
write.csv(cor_plot, "cor_plot_matrix.csv", row.names = TRUE)

# pheatmap 2
setwd("E:\\AID cohort\\code\\CellAnnotation")
library(pheatmap)
library(readxl)
df <- read_xlsx("adjusted.xlsx")
cor_adj <- as.matrix(df[, -1])
rownames(cor_adj) <- df[[1]]
#cor_adj[cor_adj < 0] <- 0
fixed_rows <- rownames(cor_adj)
fixed_cols <- colnames(cor_adj)
cor_adj <- cor_adj[fixed_rows, fixed_cols]
row_lab <- sapply(fixed_rows, function(x) strsplit(x, "_")[[1]][1])
col_lab <- fixed_cols

#color_map <- colorRampPalette(c("#0066BB", "white", "#990000"))(100)
#color_map <- colorRampPalette(c("#D2DBE8", "white", "#990000"))(100)
# color_map <- colorRampPalette(c("#FFFFE0", "#F7F7D9", "#D08886","#AE4766", "#833473", "#542D70"))(100)
color_map <- colorRampPalette(c("#FFFFE8","#FFF9D9","#FDF1C7","#FAEAC0","#F7E2BB","#F3D8B5","#EFCDAF","#EBC3A8","#E5B19E","#DEA09A","#D88F8C","#CF7D83","#C26678","#AE4766","#8E3B72","#6F326F","#542D70"))(100)


pheatmap(mat=cor_adj,color=color_map,
    cluster_rows=F,
    cluster_cols=F,
    treeheight_row=0,treeheight_col=0,
    show_rownames=T,show_colnames=T,
    cellwidth=11.8,cellheight=4.98,angle_col=45,
    labels_row=row_lab,labels_col=col_lab,
    border_color=NA,main="",legend_title="")

#convert coordinate
setwd("E:\\AID cohort\\code\\CellAnnotation")
library(pheatmap)
library(readxl)
df <- read_xlsx("adjusted.xlsx")
cor_adj <- as.matrix(df[, -1])
rownames(cor_adj) <- df[[1]]
#cor_adj[cor_adj < 0] <- 0
fixed_rows <- rownames(cor_adj)
fixed_cols <- colnames(cor_adj)
cor_adj <- t(cor_adj[fixed_rows, fixed_cols])
row_lab <- fixed_cols
col_lab <- sapply(fixed_rows, function(x) strsplit(x, "_")[[1]][1])

color_map <- colorRampPalette(c("#FFFFE8","#FFF9D9","#FDF1C7","#FAEAC0","#F7E2BB","#F3D8B5","#EFCDAF","#EBC3A8","#E5B19E","#DEA09A","#D88F8C","#CF7D83","#C26678","#AE4766","#8E3B72","#6F326F","#542D70"))(100)
pheatmap(mat=cor_adj,color=color_map,
         cluster_rows=F,
         cluster_cols=F,
         fontsize = 6,
         treeheight_row=0,treeheight_col=0,
         show_rownames=T,show_colnames=T,
         cellwidth=7.28,cellheight=14.76,angle_col=90,
         labels_row=row_lab,labels_col=col_lab,
         border_color=NA,main="",legend_title="")


#################
################# Additional Testing by GSEA analysis
library(openxlsx)
library(tidyverse)
pbmc1<-readRDS("E:\\AID cohort\\code\\result\\PBMC1.rds")
M1<-filter(pbmc1@meta.data,Level1%in%c(
  "LDG","pDC","Basophil","Platelet","Erythrocyte",
  "Naive CD4 T","Naive CD8 T"))%>%select(harmony_clusters)%>%as.matrix()%>%
  as.character()%>%unique()
Major_cell<-read.xlsx("E:\\AID cohort\\code\\result\\r5_50_t3k_AllMarkers.xlsx")
Major_cell<-Major_cell[Major_cell$cluster%in%M1,]
#unique(Major_cell$cluster)


Major_cell<-Major_cell[Major_cell$cluster%in%c(),]

Major_cell$celltype <- case_when(
  Major_cell$cluster == "0" ~ "Naive CD4 T",
  Major_cell$cluster == "4" ~ "Naive CD8 T",
  Major_cell$cluster == "8" ~ "LDG",
  Major_cell$cluster == "10" ~ "Platelet",
  Major_cell$cluster == "12" ~ "Naive CD4 T",
  Major_cell$cluster == "13" ~ "Erythrocyte",
  Major_cell$cluster == "19" ~ "Basophil",
  Major_cell$cluster == "20" ~ "pDC",
  Major_cell$cluster == "26" ~ "Erythrocyte",
  TRUE ~ "Unknown" )
unique(Major_cell$celltype)
Major_cell$cluster<-NULL



library(readxl)
library(openxlsx)
library(dplyr)
Minor_cell <- read.xlsx("E:\\AID cohort\\code\\result\\r12_AllSubCell_Combined_AllMarkers_Total.xlsx")
anno_map <- list(
  "B" = c(
    "0" = "B01_Naive_BACH2",
    "1" = "B02_Transitional_TCL1A",
    "2" = "B03_MemorySwitch_SNED1",
    "3" = "B04_CD27+BCR+",
    "4" = "B01_Naive_BACH2",
    "5" = "B05_IFNresponse_MX1",
    "6" = "B06_MemorySwitch_IGHA1"
  ),
  "DC" = c(
    "0" = "DC01_cDC2_CD1C",
    "1" = "DC02_DC_TRIM33",
    "2" = "DC03_aDC_CCR7",
    "3" = "DC04_DC_CCL5",
    "4" = "DC05_cDC2_LY96",
    "5" = "DC06_cDC1_XCR1",
    "6" = "DC07_mDC_ISG15",
    "7" = "DC08_DC_C1QTNF4",
    "8" = "DC09_mDC_VCAN"
  ),
  "NK" = c(
    "0" = "NK01_CD56dimCD16_KLRC2",
    "1" = "NK02_CD16high_PRF1",
    "2" = "NK03_CD56high_GZMK",
    "3" = "NK04_IFNresponse_MX1",
    "4" = "NK04_IFNresponse_MX1"
  ),
  "Plasma" = c(
    "0" = "Plasma01_JCHAIN",
    "1" = "Plasma02_CXCR4",
    "2" = "Plasma03_Cycling_MKI67",
    "3" = "Plasma04_XBP1"
  ),
  "NonN_CD4" = c(
    "0" = "NnCD4T01_Memory_RORA",
    "1" = "NnCD4T02_CM_FOS",
    "2" = "NnCD4T03_Treg_IKZF2",
    "3" = "NnCD4T04_CM_IFN_HERC5",
    "4" = "NnCD4T05_EMRA_GNLY",
    "5" = "NnCD4T06_CM_IFN_CISH"
  ),
  "NonN_CD8" = c(
    "0" = "NnCD8T01_Effector_ZEB2",
    "1" = "NnCD8T02_Cytotoxic_PRF1",
    "2" = "NnCD8T03_CM_NELL2",
    "3" = "NnCD8T04_EM_CXCR6",
    "4" = "NnCD8T02_Cytotoxic_PRF1",
    "5" = "NnCD8T05_EM_CMC1",
    "6" = "NnCD8T06_CM_IFN_MX1"
  ),
  "UTC" = c(
    "0" = "UTC01_MAIT_GZMK",
    "1" = "UTC01_MAIT_GZMK",
    "2" = "UTC02_gdTV2Vd9",
    "3" = "UTC03_gdTV1",
    "4" = "UTC01_MAIT_GZMK",
    "5" = "UTC03_gdTV1",
    "6" = "UTC01_MAIT_GZMK",
    "7" = "UTC01_MAIT_GZMK",
    "8" = "UTC02_gdTV2Vd9",
    "9" = "UTC01_MAIT_GZMK"
  ),
  "Mono" = c(
    "0" = "Mono01_Classical_FOSB",
    "1" = "Mono02_Inflammatory",
    "2" = "Mono03_Nonclassical",
    "3" = "Mono04_Classical_RETN",
    "4" = "Mono05_Classical_S100A9",
    "5" = "Mono06_Classical_TLN1",
    "6" = "Mono01_Classical_FOSB",
    "7" = "Mono07_IFNresponse",
    "8" = "Mono07_IFNresponse",
    "9" = "Mono02_Inflammatory"
  )
)

get_level2 <- function(ct, clu) {
  clu_str <- as.character(clu)
  map_vec <- anno_map[[ct]]
  if (!is.null(map_vec) && clu_str %in% names(map_vec)) {
    return(map_vec[[clu_str]])
  } else {
    return("Unknown_Sub")
  }
}

Minor_cell <- Minor_cell %>%
  mutate(
    Level2 = mapply(get_level2, cell_type, cluster)
  )

View(head(Minor_cell))
unique(Minor_cell$Level2)
filter(pbmc1@meta.data,Level1%in%c(
  "Monocyte","DC","UTC","B cell","Plasma","NK","Non-naive CD8 T","Non-naive CD4 T"
))%>%select(Level2)%>%as.matrix()%>%as.character()%>%unique()->m2
length(m2)
length(intersect(m2,unique(Minor_cell$Level2)))

Minor_cell$cluster<-NULL
Minor_cell$cell_type<-NULL
colnames(Minor_cell)[7]<-"celltype"
colnames(Minor_cell)
colnames(Major_cell)
AllCells_Marker<-rbind(Major_cell,Minor_cell)
write.xlsx(AllCells_Marker,file="Allcells_Marker gene.xlsx")

#GSEA analysis=====================
setwd("E:\\AID cohort\\code\\CellAnnotation")
library(openxlsx)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)

AllCells_Marker <- read.xlsx("Allcells_Marker gene.xlsx")

msig_c7 <- msigdbr(species = "Homo sapiens", category = "C7") %>%
  dplyr::select(gs_name, entrez_gene)

all_cell_types <- unique(AllCells_Marker$celltype)
gsea_all_list <- list()
gsea_pos_nes_list <- list()

for (ct in all_cell_types) {
  cat("===== processing:", ct, "=====\n")
  
  ct_markers <- AllCells_Marker %>%
    dplyr::filter(celltype == ct, avg_log2FC > 0) %>%
    dplyr::select(gene, avg_log2FC) %>%
    dplyr::arrange(dplyr::desc(avg_log2FC))
  
  if (nrow(ct_markers) < 10) {
    cat(ct, "up marker less than 10, skip\n\n")
    next
  }
  
  gene_map <- bitr(
    geneID = ct_markers$gene,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  rank_entrez <- merge(
    x = ct_markers,
    y = gene_map,
    by.x = "gene",
    by.y = "SYMBOL"
  )
  
  if (nrow(rank_entrez) < 10) {
    cat(ct, "ENTREZ matched gene less than 10, skip\n\n")
    next
  }
  
  gene_rank_vec <- rank_entrez$avg_log2FC
  names(gene_rank_vec) <- rank_entrez$ENTREZID
  gene_rank_vec <- sort(gene_rank_vec, decreasing = TRUE)
  
  gsea_res <- GSEA(
    geneList = gene_rank_vec,
    TERM2GENE = msig_c7,
    pvalueCutoff = 0.05,
    minGSSize = 10,
    maxGSSize = 5000,
    verbose = FALSE
  )
  
  gsea_all_list[[ct]] <- gsea_res
  
  gsea_clean <- gsea_res@result %>%
    dplyr::filter(NES > 0) %>%  
    dplyr::mutate(
      cell_type = ct,
      ID = sub("^GSE\\d+_", "", ID),
      Description = sub("^GSE\\d+_", "", Description),
      Rank = dplyr::row_number()
    ) %>%
    dplyr::select(cell_type, dplyr::everything())
  
  gsea_pos_nes_list[[ct]] <- gsea_clean
  
  cat(ct, "Positive enrichment pathway:", nrow(gsea_clean), "\n\n")
}


gsea_pos_all_df <- dplyr::bind_rows(gsea_pos_nes_list)#5941
gsea_pos_sig<-gsea_pos_all_df[gsea_pos_all_df$p.adjust<0.01,]

write.xlsx(
  gsea_pos_all_df,
  file = "GSEA_C7_AllCellType_PositiveNES_Result.xlsx",
  rowNames = FALSE
)
write.xlsx(
  gsea_pos_sig,
  file = "GSEA_C7_sig.xlsx",
  rowNames = FALSE
)


# if(length(gsea_pos_nes_list) > 0){
#   first_ct <- names(gsea_pos_nes_list)[1]
#   top1_pathway <- gsea_pos_nes_list[[first_ct]]$ID[1]
#   
#   gseaplot2(
#     x = gsea_all_list[[first_ct]],
#     geneSetID = top1_pathway,
#     title = paste0(first_ct, " Top1 Enrichment Pathway"),
#     color = "#2E86AB"
#   )
# }

library(readxl)
library(dplyr)
library(ggplot2)
library(scales)

library(readxl)
library(dplyr)
library(ggplot2)
library(scales)




gsea <- read_excel("GSEA_C7_sig.xlsx") %>%
  mutate(
    ID = factor(ID, levels = rev(ID)),
    logP = -log10(as.numeric(p.adjust)),
    label_x = NES + 0.15)

#barplot
#col2 <- colorRampPalette(c("#4E70AF","#7091C7","#9EBCDB","#C8D6E7","#E8EDF1","#F2EBE5","#ECD0B4","#C16D58","#A13D3B","#951d22"))(101)
col2<-c("#C0BFC9","#B6B2CC","#A89CC5","#9D8CBD",
  "#8575B1","#5E5DA0","#5A589D","#57549F","#4A539A",
  "#49549A",  "#415199") 
col2<-c("#474898",  "#5557A3","#7C76B3","#BDB8D8",  
  "#F8E7E5","#F6A5A7","#F17A80","#E55A58",  "#D33E34")


ggplot(gsea, aes(NES, ID, fill = logP)) +
  geom_col(width = 0.75) +
  geom_text(aes(x = label_x, label = cell_type), hjust = 0, size = 3.5) +
  scale_fill_gradientn(
    colours =col2,
    name="-log10(p.adjust)"
  )+
  scale_x_continuous(expand = expansion(mult = c(0,0.25))) +
  theme_classic() +
  labs(x = "NES", y = NULL) +
  theme(
    axis.ticks = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 11)
  )

#barplot

col2<- c("#F8E7E5","#F6A5A7","#E55A58","#D33E34")
col2<-c("#C0BFC9","#B6B2CC","#A89CC5","#9D8CBD",
        "#8575B1","#5E5DA0","#5A589D","#57549F","#4A539A",
        "#49549A",  "#415199") 

col2<- c("#3B7C8C","#6FA3AD","#8DB7BE","#B8D0D2","#E8EEEE",
  "#F6D7D7","#E9A6A6","#C75A62","#A93F4B"
)

ggplot(gsea, aes(x = logP, y = ID, fill = NES)) +
  geom_col(width = 0.89) +
  geom_text(aes(x = label_x, label = cell_type),
            hjust = 0, size = 3.5) +
  geom_vline(xintercept = -log10(0.05),
             linetype = "dashed",
             color = "grey60",
             linewidth = 0.6) +
  scale_fill_gradientn(
    colours =col2,
    name = "NES"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0,0.25))) +
  theme_classic() +
  labs(x = "-log10(p.adjust)", y = NULL) +
  theme(
    axis.ticks = element_line(colour="black"),
    axis.text = element_text(colour="black"),
    axis.text.y = element_text(size=10),
    axis.text.x = element_text(size=11)
  )



#lolliplot
ggplot(gsea, aes(y = ID)) +
  geom_segment(aes(x = 0, xend = logP, yend = ID),
               linewidth = 0.8, color = "grey70") +
  geom_vline(xintercept = -log10(0.01),
             linetype = "dashed",
             color = "grey60",
             linewidth = 0.6) +
  geom_point(aes(x = logP, size = NES, color = NES)) +
  scale_size_continuous(range = c(3,8), name = "NES") +
  scale_color_gradientn(
    colours = c("#474898","#BDB8D8","#F8E7E5","#E55A58","#D33E34"),
    name = "NES"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.02,0.15))) +
  theme_classic() +
  labs(x = "-log10(p.adjust)", y = NULL) +
  theme(
    axis.ticks = element_line(colour="black"),
    axis.text = element_text(colour="black"),
    axis.text.y = element_text(size=10),
    axis.text.x = element_text(size=11)
  )
