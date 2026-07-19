library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(tidyverse)
library(dplyr)

setwd("E:\\AID cohort\\code\\result")
#pbmc1<-readRDS("50_t3k_full_backup.rds")
pbmc1<-readRDS("PBMC1.rds")
Idents(pbmc1)<-"harmony_clusters"
pbmc1<-RenameIdents(pbmc1,
                    '0' = 'Naive CD4 T','1' = 'Non-naive CD8 T','2' = 'B cell',
                    '3' = 'Non-naive CD4 T','4' = 'Naive CD8 T',
                    '5' = 'Monocyte','6' = 'NK','7' = 'Monocyte','8' = 'LDG',
                    '9' = 'Monocyte',
                    '10' = 'Platelet','11' = 'Monocyte','12' = 'Naive CD4 T',
                    '13' = 'Erythrocyte','14' = 'Plasma',
                    '15' = 'Monocyte','16' = 'DC','17' = 'UTC',
                    '18' = 'Monocyte','19' = 'Basophil',
                    '20' = 'pDC','21' = 'Monocyte','22' = 'B cell',
                    '23' = 'Monocyte','24' = 'B cell',
                    '25' = 'Monocyte','26' = 'Erythrocyte','27' = 'Non-naive CD4 T',
                    '28' = 'DC','29' = 'DC',
                    '30' ='Basophil', '31' = 'Monocyte',
                    '32' ='Monocyte','33'='Monocyte')
pbmc1@meta.data$Level1 <- Idents(pbmc1)


gene_list <- list(
  "NFKB pathway" = c("TNFAIP3", "IRAK2", "NOD2", "RELA"),
  "Uncategorized" = c("ADA2", "ELF4", "CSF3R", "STAT4"),
  "Osteoclast function" = c("ACP5", "OGFRL1"),
  "Ca2+_flux-PLC" = c("PLCG2", "PLCG1"),
  "Endolysosomal_nucleic_acid_sensing" = c("UNC93B1", "TLR8", "TLR7", "TLR1", "PLD4", "TRAF3"),
  "Arachidonic acid metabolism" = c("TBXAS1"),
  "Inflammasome_IL-1β" = c("IL1R1", "NLRC4", "NLRP3", "PSTPIP1", "CDC42", "LPIN2"),
  "Immune_metabolic" = c("LACC1", "SLC7A7"),
  "PP" = c("polygenic"),
  "HC" = c("wide-type"),
  "Inborn_errors_of_cell_death" = c("RIPK1", "TNFRSF1A", "OTULIN", "RNF31"),
  "Negative_regulation_of_IFN-I" = c("ISG15", "USP18", "SOCS1"),
  "Cytoskeleton_and_small_GTPase" = c("GNAI2", "KRAS"),
  "TBK1_IRF3" = c("STING", "TREX1", "IFIH1", "COPA"),
  "Protein_homeostasis" = c("UBA1", "PSMD12"))

#####1 Ro/e calculation (heatmap and scatter plots)
library(dplyr)
library(stringr)
library(tidyverse)

####(1) heatmap
plot_RoE_heatmap <- function(seurat_obj, gene_list, group_type = c("condition", "gene")){
  group_type <- match.arg(group_type)
  meta <- seurat_obj@meta.data
  
  if(group_type == "condition"){
    gene_map <- stack(gene_list)
    colnames(gene_map) <- c("gene", "condition0")
    meta$condition0 <- gene_map$condition0[match(meta$gene, gene_map$gene)]
    
    meta <- meta %>% mutate(
      condition = case_when(
        str_detect(samples, "SLE") ~ "SLE",
        str_detect(samples, "JIA") ~ "JIA",
        str_detect(samples, "KD") ~ "KD",
        str_detect(samples, "HC") ~ "HC",
        TRUE ~ condition0
      )
    )
    tbl <- table(meta$condition, meta$Level1)
    row_var <- "condition"
    target_last <- c("SLE", "JIA", "KD", "HC")
  }else{
    gene_map <- stack(gene_list)
    colnames(gene_map) <- c("gene", "pathway")
    meta$pathway <- gene_map$pathway[match(meta$gene, gene_map$gene)]
    
    meta <- meta %>% mutate(gene2 = gene)
    meta <- meta %>% mutate(
      gene2 = case_when(
        str_detect(samples, "SLE") ~ "SLE",
        str_detect(samples, "JIA") ~ "JIA",
        str_detect(samples, "KD") ~ "KD",
        str_detect(samples, "HC") ~ "HC",
        TRUE ~ gene2
      )
    )
    tbl <- table(meta$gene2, meta$Level1)
    row_var <- "gene2"
    target_last <- c("SLE", "JIA", "KD", "HC")
  }
  
  obs_mat <- as.data.frame.matrix(tbl)
  total <- sum(obs_mat)
  row_sum <- rowSums(obs_mat)
  col_sum <- colSums(obs_mat)
  E <- outer(row_sum, col_sum) / total
  RoE <- obs_mat / E
  
  cell_order <- names(sort(table(meta$Level1), decreasing = TRUE))
  RoE_fix <- RoE[, cell_order]
  
  grade_mat <- matrix(nrow = nrow(RoE_fix), ncol = ncol(RoE_fix))
  rownames(grade_mat) <- rownames(RoE_fix)
  colnames(grade_mat) <- colnames(RoE_fix)
  for(i in seq_len(nrow(RoE_fix))){
    for(j in seq_len(ncol(RoE_fix))){
      val <- RoE_fix[i,j]
      grade_mat[i,j] <- case_when(
        val <= 1 ~ "-",
        val > 1 & val <= 2 ~ "+",
        val > 2 & val <= 4 ~ "++",
        val > 4 ~ "+++"
      )
    }
  }
  
  if(group_type == "condition"){
    all_rows <- setdiff(rownames(RoE_fix), target_last)
    row_order <- c(all_rows, target_last)
  }else{
    gene_path_df <- distinct(meta, gene2, pathway) %>% filter(!gene2 %in% target_last)
    sorted_gene <- gene_path_df %>% arrange(pathway, gene2) %>% pull(gene2)
    row_order <- c(sorted_gene, target_last)
  }
  RoE_final <- RoE_fix[row_order, ]
  grade_final <- grade_mat[row_order, ]
  
  df_roe <- as.data.frame(RoE_final) %>%
    rownames_to_column(row_var) %>%
    pivot_longer(-all_of(row_var), names_to = "celltype", values_to = "RoE")
  df_label <- as.data.frame(grade_final) %>%
    rownames_to_column(row_var) %>%
    pivot_longer(-all_of(row_var), names_to = "celltype", values_to = "label")
  plot_df <- inner_join(df_roe, df_label, by = c(row_var, "celltype"))
  
  plot_df <- plot_df %>% mutate(
    roe_level = case_when(
      RoE <= 1 ~ "0-1",
      RoE > 1 & RoE <= 2 ~ "1-2",
      RoE > 2 & RoE <= 4 ~ "2-4",
      RoE > 4 ~ ">4"
    ),
    roe_level = factor(roe_level, levels = c(">4", "2-4", "1-2", "0-1"))
  )
  
  level_color <- c(
    ">4" = "#912623",
    "2-4" = "#E55C3A",
    "1-2" = "#F7BE6F",
    "0-1" = "#FFEFC1"
  )
  
  plot_df$celltype <- factor(plot_df$celltype, levels = cell_order)
  plot_df[[row_var]] <- factor(plot_df[[row_var]], levels = rev(row_order))
  
  p <- ggplot(plot_df, aes(x = celltype, y = .data[[row_var]])) +
    geom_tile(aes(fill = roe_level), width = 0.98, height = 0.98) +
    geom_text(aes(label = label), colour = "black", size = 3) +
    scale_fill_manual(values = level_color) +
    labs(fill = "Ro/e", x = NULL, y = NULL) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
      axis.text.y = element_text(colour = "black"),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      legend.text = element_text(colour = "black"),
      legend.title = element_text(colour = "black")
    )
  
  res_list <- list(
    figure = p,
    RoE_matrix = RoE_final
  )
  return(res_list)
}

# Test
out_cond <- plot_RoE_heatmap(pbmc1, gene_list, group_type = "condition")
out_cond$figure
out_gene <- plot_RoE_heatmap(pbmc1, gene_list, group_type = "gene")
out_gene$figure


####(2) scatter
plot_RoE_scatter <- function(seurat_obj, gene_list, 
group_type = c("condition", "gene"), 
point_size = 2, legend_key_size = 0.5){
  group_type <- match.arg(group_type)
  meta <- seurat_obj@meta.data
  
  my_colors <- c(
    "Naive CD4 T"="#847AB3","Non-naive CD4 T"="#C8BBD5",
    "Naive CD8 T"="#F6E1EE",
    "Non-naive CD8 T"="#0B71AB","NK"="#96AF95",
    "B cell"="#003D81","Monocyte"="#F49D5C",
    "LDG"="#C95968","Plasma"="#B3D1E7",
    "Erythrocyte"="#FFEFC1","Platelet"="#B3928B",
    "UTC"="#86C7B4","DC"="#F0E442",
    "pDC"="#B07A99","Basophil"="#8C2522"
  )
  
  if(group_type == "condition"){
    gene_map <- stack(gene_list)
    colnames(gene_map) <- c("gene", "condition0")
    meta$condition0 <- gene_map$condition0[match(meta$gene, gene_map$gene)]
    
    meta <- meta %>% mutate(
      condition = case_when(
        str_detect(samples, "SLE") ~ "SLE",
        str_detect(samples, "JIA") ~ "JIA",
        str_detect(samples, "KD") ~ "KD",
        str_detect(samples, "HC") ~ "HC",
        TRUE ~ condition0
      )
    )
    tbl <- table(meta$condition, meta$Level1)
    row_var <- "condition"
    target_last <- c("SLE", "JIA", "KD", "HC")
  }else{
    gene_map <- stack(gene_list)
    colnames(gene_map) <- c("gene", "pathway")
    meta$pathway <- gene_map$pathway[match(meta$gene, gene_map$gene)]
    
    meta <- meta %>% mutate(gene2 = gene)
    meta <- meta %>% mutate(
      gene2 = case_when(
        str_detect(samples, "SLE") ~ "SLE",
        str_detect(samples, "JIA") ~ "JIA",
        str_detect(samples, "KD") ~ "KD",
        str_detect(samples, "HC") ~ "HC",
        TRUE ~ gene2
      )
    )
    tbl <- table(meta$gene2, meta$Level1)
    row_var <- "gene2"
    target_last <- c("SLE", "JIA", "KD", "HC")
  }
  
  obs_mat <- as.data.frame.matrix(tbl)
  total <- sum(obs_mat)
  row_sum <- rowSums(obs_mat)
  col_sum <- colSums(obs_mat)
  E <- outer(row_sum, col_sum) / total
  RoE <- obs_mat / E
  
  cell_order <- names(sort(table(meta$Level1), decreasing = TRUE))
  RoE_fix <- RoE[, cell_order]
  
  if(group_type == "condition"){
    all_rows <- setdiff(rownames(RoE_fix), target_last)
    row_order <- c(all_rows, target_last)
  }else{
    gene_path_df <- distinct(meta, gene2, pathway) %>% filter(!gene2 %in% target_last)
    sorted_gene <- gene_path_df %>% arrange(pathway, gene2) %>% pull(gene2)
    row_order <- c(sorted_gene, target_last)
  }
  RoE_final <- RoE_fix[row_order, ]
  
  df_roe <- as.data.frame(RoE_final) %>%
    rownames_to_column(row_var) %>%
    pivot_longer(-all_of(row_var), names_to = "celltype", values_to = "RoE")
  
  df_roe$celltype <- factor(df_roe$celltype, levels = cell_order)
  df_roe[[row_var]] <- factor(df_roe[[row_var]], levels = row_order)
  
  if(group_type == "condition"){
    x_text_theme <- element_text(angle = 45, hjust = 1, colour = "black")
  }else{
    x_text_theme <- element_text(angle = 90, hjust = 1, colour = "black")
  }
  
  p <- ggplot(df_roe, aes(x = .data[[row_var]], y = RoE, colour = celltype)) +
    geom_point(size = point_size) +
    geom_hline(yintercept = c(0,1,2,4,6), linetype = "solid", linewidth = 0.3) +
    scale_y_continuous(breaks = c(0,1,2,4,6), limits = c(0,6)) +
    scale_colour_manual(values = my_colors, drop = FALSE) +
    labs(x = NULL, y = "Ro/e", colour = "Cell type") +
    theme_classic() +
    theme(
      axis.text.x = x_text_theme,
      axis.text.y = element_text(colour = "black"),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      legend.text = element_text(colour = "black"),
      legend.title = element_text(colour = "black")
    ) +
    guides(colour = guide_legend(reverse = FALSE, keywidth = legend_key_size, keyheight = legend_key_size, default.unit = "cm"))
  return(p)
}
plot_RoE_scatter(pbmc1, gene_list, group_type = "condition")

######2 Prop calculation
library(dplyr)
library(stringr)
library(tidyverse)

####(1) stack
plot_stack_bar_freq <- function(seurat_obj, gene_list, group_type = c("condition", "gene"), 
                          legend_key_size = 0.5){
  group_type <- match.arg(group_type)
  meta <- seurat_obj@meta.data
  
  my_colors <- c(
    "Naive CD4 T"="#847AB3","Non-naive CD4 T"="#C8BBD5",
    "Naive CD8 T"="#F6E1EE",
    "Non-naive CD8 T"="#0B71AB","NK"="#96AF95",
    "B cell"="#003D81","Monocyte"="#F49D5C",
    "LDG"="#C95968","Plasma"="#B3D1E7",
    "Erythrocyte"="#FFEFC1","Platelet"="#B3928B",
    "UTC"="#86C7B4","DC"="#F0E442",
    "pDC"="#B07A99","Basophil"="#8C2522"
  )
  
  if(group_type == "condition"){
    gene_map <- stack(gene_list)
    colnames(gene_map) <- c("gene", "condition0")
    meta$condition0 <- gene_map$condition0[match(meta$gene, gene_map$gene)]
    
    meta <- meta %>% mutate(
      condition = case_when(
        str_detect(samples, "SLE") ~ "SLE",
        str_detect(samples, "JIA") ~ "JIA",
        str_detect(samples, "KD") ~ "KD",
        str_detect(samples, "HC") ~ "HC",
        TRUE ~ condition0
      )
    )
    tbl <- table(meta$condition, meta$Level1)
    row_var <- "condition"
    target_last <- c("SLE", "JIA", "KD", "HC")
  }else{
    gene_map <- stack(gene_list)
    colnames(gene_map) <- c("gene", "pathway")
    meta$pathway <- gene_map$pathway[match(meta$gene, gene_map$gene)]
    
    meta <- meta %>% mutate(gene2 = gene)
    meta <- meta %>% mutate(
      gene2 = case_when(
        str_detect(samples, "SLE") ~ "SLE",
        str_detect(samples, "JIA") ~ "JIA",
        str_detect(samples, "KD") ~ "KD",
        str_detect(samples, "HC") ~ "HC",
        TRUE ~ gene2
      )
    )
    tbl <- table(meta$gene2, meta$Level1)
    row_var <- "gene2"
    target_last <- c("SLE", "JIA", "KD", "HC")
  }
  
  obs_mat <- as.data.frame.matrix(tbl)
  cell_order <- names(sort(table(meta$Level1), decreasing = TRUE))
  obs_mat <- obs_mat[, cell_order]
  
  Freq_mat <- as.data.frame(t(apply(obs_mat,1,function(x){
    y <- as.numeric(x)
    y1 <- y / sum(y)
    return(y1)
  })))
  colnames(Freq_mat) <- colnames(obs_mat)
  
  if(group_type == "condition"){
    all_rows <- setdiff(rownames(Freq_mat), target_last)
    row_order <- c(all_rows, target_last)
  }else{
    gene_path_df <- distinct(meta, gene2, pathway) %>% filter(!gene2 %in% target_last)
    sorted_gene <- gene_path_df %>% arrange(pathway, gene2) %>% pull(gene2)
    row_order <- c(sorted_gene, target_last)
  }
  Freq_final <- Freq_mat[row_order, ]
  
  plot_df <- as.data.frame(Freq_final) %>%
    rownames_to_column(row_var) %>%
    pivot_longer(-all_of(row_var), names_to = "celltype", values_to = "Freq")
  
  plot_df$celltype <- factor(plot_df$celltype, levels = cell_order)
  plot_df[[row_var]] <- factor(plot_df[[row_var]], levels = row_order)
  
  if(group_type == "condition"){
    x_text_theme <- element_text(angle = 90, hjust = 1, colour = "black")
  }else{
    x_text_theme <- element_text(angle = 90, hjust = 1, colour = "black")
  }
  
  p <- ggplot(plot_df, aes(x = .data[[row_var]], y = Freq, fill = celltype)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = my_colors, drop = FALSE) +
    labs(x = NULL, y = "Cell proportion", fill = "Cell type") +
    theme_classic() +
    theme(
      axis.text.x = x_text_theme,
      axis.text.y = element_text(colour = "black"),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      legend.text = element_text(colour = "black"),
      legend.title = element_text(colour = "black")
    ) +
    guides(fill = guide_legend(reverse = FALSE, keywidth = legend_key_size, keyheight = legend_key_size, default.unit = "cm"))
  return(p)
}

# test run
plot_stack_bar_freq(pbmc1, gene_list, group_type = "condition", legend_key_size = 0.3)
plot_stack_bar_freq(pbmc1, gene_list, group_type = "gene", legend_key_size = 0.3)

####(2) scatter
plot_prop_scatter <- function(seurat_obj, gene_list, 
group_type = c("condition", "gene"),point_size = 2,
legend_key_size = 0.5){
  group_type <- match.arg(group_type)
  meta <- seurat_obj@meta.data
  
  my_colors <- c(
    "Naive CD4 T"="#847AB3","Non-naive CD4 T"="#C8BBD5",
    "Naive CD8 T"="#F6E1EE",
    "Non-naive CD8 T"="#0B71AB","NK"="#96AF95",
    "B cell"="#003D81","Monocyte"="#F49D5C",
    "LDG"="#C95968","Plasma"="#B3D1E7",
    "Erythrocyte"="#FFEFC1","Platelet"="#B3928B",
    "UTC"="#86C7B4","DC"="#F0E442",
    "pDC"="#B07A99","Basophil"="#8C2522"
  )
  
  if(group_type == "condition"){
    gene_map <- stack(gene_list)
    colnames(gene_map) <- c("gene", "condition0")
    meta$condition0 <- gene_map$condition0[match(meta$gene, gene_map$gene)]
    
    meta <- meta %>% mutate(
      condition = case_when(
        str_detect(samples, "SLE") ~ "SLE",
        str_detect(samples, "JIA") ~ "JIA",
        str_detect(samples, "KD") ~ "KD",
        str_detect(samples, "HC") ~ "HC",
        TRUE ~ condition0
      )
    )
    tbl <- table(meta$condition, meta$Level1)
    row_var <- "condition"
    target_last <- c("SLE", "JIA", "KD", "HC")
  }else{
    gene_map <- stack(gene_list)
    colnames(gene_map) <- c("gene", "pathway")
    meta$pathway <- gene_map$pathway[match(meta$gene, gene_map$gene)]
    
    meta <- meta %>% mutate(gene2 = gene)
    meta <- meta %>% mutate(
      gene2 = case_when(
        str_detect(samples, "SLE") ~ "SLE",
        str_detect(samples, "JIA") ~ "JIA",
        str_detect(samples, "KD") ~ "KD",
        str_detect(samples, "HC") ~ "HC",
        TRUE ~ gene2
      )
    )
    tbl <- table(meta$gene2, meta$Level1)
    row_var <- "gene2"
    target_last <- c("SLE", "JIA", "KD", "HC")
  }
  
  obs_mat <- as.data.frame.matrix(tbl)
  cell_order <- names(sort(table(meta$Level1), decreasing = TRUE))
  obs_mat <- obs_mat[, cell_order]
  
  Freq_mat <- as.data.frame(t(apply(obs_mat,1,function(x){
    y <- as.numeric(x)
    y1 <- y / sum(y)
    return(y1)
  })))
  colnames(Freq_mat) <- colnames(obs_mat)
  
  if(group_type == "condition"){
    all_rows <- setdiff(rownames(Freq_mat), target_last)
    row_order <- c(all_rows, target_last)
  }else{
    gene_path_df <- distinct(meta, gene2, pathway) %>% filter(!gene2 %in% target_last)
    sorted_gene <- gene_path_df %>% arrange(pathway, gene2) %>% pull(gene2)
    row_order <- c(sorted_gene, target_last)
  }
  Freq_final <- Freq_mat[row_order, ]
  
  plot_df <- as.data.frame(Freq_final) %>%
    rownames_to_column(row_var) %>%
    pivot_longer(-all_of(row_var), names_to = "celltype", values_to = "Freq")
  
  plot_df$celltype <- factor(plot_df$celltype, levels = cell_order)
  plot_df[[row_var]] <- factor(plot_df[[row_var]], levels = row_order)
  
  if(group_type == "condition"){
    x_text_theme <- element_text(angle = 45, hjust = 1, colour = "black")
  }else{
    x_text_theme <- element_text(angle = 90, hjust = 1, colour = "black")
  }
  
  p <- ggplot(plot_df, aes(x = .data[[row_var]], y = Freq, colour = celltype)) +
    geom_point(size = point_size) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
    scale_colour_manual(values = my_colors, drop = FALSE) +
    labs(x = NULL, y = "Cell proportion", colour = "Cell type") +
    theme_classic() +
    theme(
      axis.text.x = x_text_theme,
      axis.text.y = element_text(colour = "black"),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      legend.text = element_text(colour = "black"),
      legend.title = element_text(colour = "black")
    ) +
    guides(colour = guide_legend(reverse = FALSE, keywidth = legend_key_size, keyheight = legend_key_size, default.unit = "cm"))
  return(p)
}
plot_prop_scatter(pbmc1, gene_list, group_type = "condition", 
                legend_key_size = 0.3)



#####3 Roe and Prop merge to a point plots (Level1 and Level2)
library(dplyr)
library(stringr)
library(tidyverse)
library(scales)
#color_palette <- colorRampPalette(c("#E6F7E6", "#B2E8D8", "#80D8D8", "#3388BB", "#0055BB"))(20)
#color_palette2 <- colorRampPalette(c("#F0F9E8","#D9F0D3","#BDE0C0","#9FD1B0","#7FC2B0","#60B3B3","#40A3B6","#2093B9","#0083BC","#0070AA","#005E99","#004C88"))(100)

pal_orrd <- brewer.pal(n = 9, name = "OrRd")
pal_orrd_100 <- colorRampPalette(pal_orrd)(100)

rdbu_brewer <- brewer.pal(11,"RdBu")
rdbu_100 <- colorRampPalette(rdbu_brewer)(100)
#zscore_pal <- colorRampPalette(c("#2066BB","#70B2DD","#F8F8F8","#F7A068","#B92222"))(100)

plot_RoE_prop_square_bubble <- function(seurat_obj,
                                        gene_list,
                                        group_type = c("condition", "gene"),
                                        base_point_size = 6,
                                        legend_key_size = 0.5,
                                        color_vec = NULL) {
  group_type <- match.arg(group_type)
  meta <- seurat_obj@meta.data
  
  cell_custom_order <- c(
    "Monocyte", "LDG", "pDC", "DC", "Basophil",
    "Naive CD4 T", "Non-naive CD4 T", "Naive CD8 T", "Non-naive CD8 T", "UTC",
    "B cell", "Plasma",
    "NK",
    "Platelet", "Erythrocyte"
  )
  target_last <- c("SLE", "JIA", "KD", "HC")
  fixed_pathway_order <- c(
    "NFKB pathway",
    "Osteoclast function",
    "Ca2+_flux-PLC",
    "Endolysosomal_nucleic_acid_sensing",
    "Arachidonic acid metabolism",
    "Inflammasome_IL-1β",
    "Immune_metabolic",
    "Inborn_errors_of_cell_death",
    "Negative_regulation_of_IFN-I",
    "Cytoskeleton_and_small_GTPase",
    "TBK1_IRF3",
    "Protein_homeostasis",
    "Uncategorized"
  )
  
  if (is.null(color_vec)) {
    prop_color <- colorRampPalette(c("#E6F7E6", "#B2E8D8", "#80D8D8", "#3388BB", "#0055BB"))(100)
  } else {
    prop_color <- colorRampPalette(color_vec)(100)
  }
  
  if(group_type == "condition"){
    gene_map <- stack(gene_list)
    colnames(gene_map) <- c("gene", "condition0")
    meta$condition0 <- gene_map$condition0[match(meta$gene, gene_map$gene)]
    
    meta <- meta %>% mutate(
      condition = case_when(
        str_detect(samples, "SLE") ~ "SLE",
        str_detect(samples, "JIA") ~ "JIA",
        str_detect(samples, "KD") ~ "KD",
        str_detect(samples, "HC") ~ "HC",
        TRUE ~ condition0
      )
    )
    tbl <- table(meta$condition, meta$Level1)
    row_var <- "condition"
  }else{
    gene_map <- stack(gene_list)
    colnames(gene_map) <- c("gene", "pathway")
    meta$pathway <- gene_map$pathway[match(meta$gene, gene_map$gene)]
    
    meta <- meta %>% mutate(gene2 = gene)
    meta <- meta %>% mutate(
      gene2 = case_when(
        str_detect(samples, "SLE") ~ "SLE",
        str_detect(samples, "JIA") ~ "JIA",
        str_detect(samples, "KD") ~ "KD",
        str_detect(samples, "HC") ~ "HC",
        TRUE ~ gene2
      )
    )
    tbl <- table(meta$gene2, meta$Level1)
    row_var <- "gene2"
  }
  
  obs_mat <- as.data.frame.matrix(tbl)
  total <- sum(obs_mat)
  row_sum <- rowSums(obs_mat)
  col_sum <- colSums(obs_mat)
  E <- outer(row_sum, col_sum) / total
  RoE <- obs_mat / E
  
  cell_order <- intersect(cell_custom_order, colnames(obs_mat))
  RoE_fix <- RoE[, cell_order]
  
  Freq_mat <- as.data.frame(t(apply(obs_mat, 1, function(x) {
    y <- as.numeric(x)
    y1 <- y / sum(y)
    return(y1)
  })))
  colnames(Freq_mat) <- colnames(obs_mat)
  Freq_fix <- Freq_mat[, cell_order]
  
  grade_mat <- matrix(nrow = nrow(RoE_fix), ncol = ncol(RoE_fix))
  rownames(grade_mat) <- rownames(RoE_fix)
  colnames(grade_mat) <- colnames(RoE_fix)
  for(i in seq_len(nrow(RoE_fix))){
    for(j in seq_len(ncol(RoE_fix))){
      val <- RoE_fix[i,j]
      grade_mat[i,j] <- case_when(
        val <= 1 ~ "-",
        val > 1 & val <= 2 ~ "+",
        val > 2 & val <= 4 ~ "++",
        val > 4 ~ "+++"
      )
    }
  }
  
  all_rows <- rownames(RoE_fix)
  disease_rows <- intersect(target_last, all_rows)
  if(group_type == "condition"){
    pathway_rows <- intersect(fixed_pathway_order, all_rows)
    row_order <- c(pathway_rows, disease_rows)
  }else{
    gene_path_df <- distinct(meta, gene2, pathway) %>% filter(!gene2 %in% target_last)
    gene_path_df$pathway <- factor(gene_path_df$pathway, levels = fixed_pathway_order)
    sorted_gene <- gene_path_df %>% arrange(pathway, gene2) %>% pull(gene2)
    row_order <- c(sorted_gene, disease_rows)
  }
  
  RoE_final <- RoE_fix[row_order, ]
  Freq_final <- Freq_fix[row_order, ]
  grade_final <- grade_mat[row_order, ]
  
  df_roe <- as.data.frame(RoE_final) %>%
    rownames_to_column(row_var) %>%
    pivot_longer(-all_of(row_var), names_to = "celltype", values_to = "RoE")
  df_freq <- as.data.frame(Freq_final) %>%
    rownames_to_column(row_var) %>%
    pivot_longer(-all_of(row_var), names_to = "celltype", values_to = "Freq")
  
  plot_df <- df_roe %>%
    inner_join(df_freq, by = c(row_var, "celltype")) %>%
    mutate(
      roe_level = case_when(
        RoE <= 1 ~ "0-1",
        RoE > 1 & RoE <= 2 ~ "1-2",
        RoE > 2 & RoE <= 4 ~ "2-4",
        RoE > 4 ~ ">4"
      ),
      roe_level = factor(roe_level, levels = c("0-1", "1-2", "2-4", ">4")),
      size_grade = case_when(
        roe_level == "0-1"  ~ base_point_size * 0.3,
        roe_level == "1-2"  ~ base_point_size * 0.6,
        roe_level == "2-4"  ~ base_point_size * 0.9,
        roe_level == ">4"   ~ base_point_size * 1.2
      )
    )
  
  plot_df$celltype <- factor(plot_df$celltype, levels = cell_order)
  plot_df[[row_var]] <- factor(plot_df[[row_var]], levels = row_order)
  
  x_text_theme <- element_text(angle = 45, hjust = 1, colour = "black")
  
  p <- ggplot(plot_df, aes(x = celltype, y = .data[[row_var]])) +
    geom_point(aes(size = size_grade, fill = Freq), shape = 22, stroke = 0) +
    scale_fill_gradientn(colours = prop_color) +
    scale_size(
      name = "Ro/e level",
      breaks = c(base_point_size*0.3, base_point_size*0.6, base_point_size*0.9, base_point_size*1.2),
      labels = c("0-1", "1-2", "2-4", ">4"),
      range = range(plot_df$size_grade)
    ) +
    labs(fill = "Cell proportion", x = NULL, y = NULL) +
    guides(
      fill = guide_colorbar(barwidth = legend_key_size, barheight = legend_key_size * 5),
      size = guide_legend(keywidth = legend_key_size, keyheight = legend_key_size, reverse = TRUE)
    ) +
    theme_classic() +
    theme(
      axis.text.x = x_text_theme,
      axis.text.y = element_text(colour = "black"),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      legend.text = element_text(colour = "black"),
      legend.title = element_text(colour = "black")
    )
  
  res_list <- list(
    figure = p,
    RoE_matrix = RoE_final,
    Prop_matrix = Freq_final
  )
  return(res_list)
}

out_condition <- plot_RoE_prop_square_bubble(
  seurat_obj = pbmc1,
  gene_list = gene_list,
  group_type = "condition",
  base_point_size = 6,
  legend_key_size = 0.5,
  color_vec = pal_orrd_100
)
out_condition$figure

out_gene <- plot_RoE_prop_square_bubble(
  seurat_obj = pbmc1,
  gene_list = gene_list,
  group_type = "gene",
  base_point_size = 5,
  legend_key_size = 0.5,
  color_vec = pal_orrd_100
)
out_gene$figure

#####4 boxplots: proportion of major celltyps differed by sample status [Level 1]
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggsignif)


status_order <- c("unaffected control", "polygenic patient", "monogenic patient")
status_col <- c("unaffected control"="#4375A0","polygenic patient"="#B18749","monogenic patient"="#912623")
cell_order <- names(sort(table(pbmc1$Level1),decreasing = T))
dodge_pos <- position_dodge(width = 0.6)

meta_df <- pbmc1@meta.data%>%
  group_by(samples) %>%
  mutate(total_sample_cells = n()) %>%
  ungroup() %>%
  count(samples, status, Level1, total_sample_cells) %>%
  mutate(cell_prop = n / total_sample_cells) %>%
  select(samples, status, Level1, cell_prop) %>%
  mutate(status = factor(status, levels = status_order),
  Level1 = factor(Level1, levels = cell_order))


kw_res <- meta_df %>%
  group_by(Level1) %>%
  group_modify(~tibble(kw_p = kruskal.test(cell_prop ~ status, data = .x)$p.value)) %>%
  ungroup() %>%
  mutate(signif_label = case_when(
    kw_p < 0.001 ~ "***",
    kw_p < 0.01 ~ "**",
    kw_p < 0.05 ~ "*",
    TRUE ~ "ns"
  ))%>%as.data.frame()%>%filter(kw_p<0.05)
kw_res<-kw_res[-grep("Erythrocyte|DC|Baso",kw_res$Level1),]
meta_df2<-meta_df[meta_df$Level1%in%kw_res$Level1,]


prop_statusbox <- ggplot(meta_df2, aes(x = Level1, y = cell_prop, color = status)) +
  geom_boxplot(position = dodge_pos, fill = NA, width = 0.55, 
               linewidth = 0.6, outlier.shape = NA) +
  scale_color_manual(values = status_col) +
  labs(x = "", y = "Cell proportion per sample") +
  theme_classic() +
  theme(
    legend.position="none",
    axis.text.x = element_text(angle = 35, hjust = 1, colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line("black"),
    panel.border = element_blank(),
    axis.line = element_line("black"),
    legend.title = element_blank(),
    legend.text = element_text(colour = "black")
  )
prop_statusbox 
options(scipen = 999)
as.numeric(kw_res$kw_p)


##############Level2  
##################################
##################################
######1  just change Level1 to Level2
library(dplyr)
library(stringr)
library(tidyverse)
library(scales)
library(RColorBrewer)

pal_orrd <- brewer.pal(n = 9, name = "OrRd")
pal_orrd_100 <- colorRampPalette(pal_orrd)(100)
rdbu_brewer <- brewer.pal(11,"RdBu")
rdbu_100 <- colorRampPalette(rdbu_brewer)(100)

plot_RoE_prop_square_bubble <- function(seurat_obj,
                                        gene_list,
                                        group_type = c("condition", "gene"),
                                        freq_mode = c("count", "norm"),
                                        base_point_size = 6,
                                        legend_key_size = 0.5,
                                        color_vec = NULL) {
  group_type <- match.arg(group_type)
  freq_mode <- match.arg(freq_mode)
  meta <- seurat_obj@meta.data
  
  cell_custom_order <- c(
    "DC01_cDC2_CD1C","DC02_DC_TRIM33","DC03_aDC_CCR7","DC04_DC_CCL5","DC05_cDC2_LY96","DC06_cDC1_XCR1","DC07_mDC_ISG15","DC08_DC_C1QTNF4","DC09_mDC_VCAN",
    "Mono01_Classical_FOSB","Mono02_Inflammatory","Mono03_Nonclassical","Mono04_Classical_RETN","Mono05_Classical_S100A9","Mono06_Classical_TLN1","Mono07_IFNresponse",
    "pDC",
    "LDG",
    "Basophil",
    "NnCD4T01_Memory_RORA","NnCD4T02_CM_FOS","NnCD4T03_Treg_IKZF2","NnCD4T04_CM_IFN_HERC5","NnCD4T05_EMRA_GNLY","NnCD4T06_CM_IFN_CISH",
    "NnCD8T01_Effector_ZEB2","NnCD8T02_Cytotoxic_PRF1","NnCD8T03_CM_NELL2","NnCD8T04_EM_CXCR6","NnCD8T05_EM_CMC1","NnCD8T06_CM_IFN_MX1",
    "B01_Naive_BACH2","B02_Transitional_TCL1A","B03_MemorySwitch_SNED1","B04_CD27+BCR+","B05_IFNresponse_MX1","B06_MemorySwitch_IGHA1",
    "Plasma01_JCHAIN","Plasma02_CXCR4","Plasma03_Cycling_MKI67","Plasma04_XBP1",
    "NK01_CD56dimCD16_KLRC2","NK02_CD16high_PRF1","NK03_CD56high_GZMK","NK04_IFNresponse_MX1",
    "UTC01_MAIT_GZMK","UTC02_gdTV2Vd9","UTC03_gdTV1",
    "Erythrocyte",
    "Platelet"
  )
  
  target_last <- c("SLE", "JIA", "KD", "HC")
  fixed_pathway_order <- c(
    "NFKB pathway",
    "Osteoclast function",
    "Ca2+_flux-PLC",
    "Endolysosomal_nucleic_acid_sensing",
    "Arachidonic acid metabolism",
    "Inflammasome_IL-1β",
    "Immune_metabolic",
    "Inborn_errors_of_cell_death",
    "Negative_regulation_of_IFN-I",
    "Cytoskeleton_and_small_GTPase",
    "TBK1_IRF3",
    "Protein_homeostasis",
    "Uncategorized"
  )
  
  if (is.null(color_vec)) {
    prop_color <- colorRampPalette(c("#E6F7E6", "#B2E8D8", "#80D8D8", "#3388BB", "#0055BB"))(100)
  } else {
    prop_color <- colorRampPalette(color_vec)(100)
  }
  
  if (group_type == "condition") {
    gene_map <- stack(gene_list)
    colnames(gene_map) <- c("gene", "condition0")
    meta$condition0 <- gene_map$condition0[match(meta$gene, gene_map$gene)]
    
    meta <- meta %>% mutate(
      condition = case_when(
        str_detect(samples, "SLE") ~ "SLE",
        str_detect(samples, "JIA") ~ "JIA",
        str_detect(samples, "KD") ~ "KD",
        str_detect(samples, "HC") ~ "HC",
        TRUE ~ condition0
      )
    )
    tbl <- table(meta$condition, meta$Level2) 
    row_var <- "condition"
  } else {
    gene_map <- stack(gene_list)
    colnames(gene_map) <- c("gene", "pathway")
    meta$pathway <- gene_map$pathway[match(meta$gene, gene_map$gene)]
    
    meta <- meta %>% mutate(gene2 = gene)
    meta <- meta %>% mutate(
      gene2 = case_when(
        str_detect(samples, "SLE") ~ "SLE",
        str_detect(samples, "JIA") ~ "JIA",
        str_detect(samples, "KD") ~ "KD",
        str_detect(samples, "HC") ~ "HC",
        TRUE ~ gene2
      )
    )
    tbl <- table(meta$gene2, meta$Level2) 
    row_var <- "gene2"
  }
  
  obs_mat <- as.data.frame.matrix(tbl)
  total <- sum(obs_mat)
  row_sum <- rowSums(obs_mat)
  col_sum <- colSums(obs_mat)
  E <- outer(row_sum, col_sum) / total
  RoE <- obs_mat / E
  
  cell_order <- intersect(cell_custom_order, colnames(obs_mat))
  RoE_fix <- RoE[, cell_order]
  

  Freq_raw_mat <- as.data.frame(t(apply(obs_mat, 1, function(x) {
    y <- as.numeric(x)
    y1 <- y / sum(y)
    return(y1)
  })))
  colnames(Freq_raw_mat) <- colnames(obs_mat)
  Freq_raw_fix <- Freq_raw_mat[, cell_order]
  
  Freq_norm_fix <- Freq_raw_fix
  for (cell in colnames(Freq_norm_fix)) {
    vec <- Freq_raw_fix[[cell]]
    minv <- min(vec, na.rm = T)
    maxv <- max(vec, na.rm = T)
    if (maxv - minv > 1e-8) {
      Freq_norm_fix[[cell]] <- (vec - minv) / (maxv - minv)
    } else {
      Freq_norm_fix[[cell]] <- 0
    }
  }
  
  grade_mat <- matrix(nrow = nrow(RoE_fix), ncol = ncol(RoE_fix))
  rownames(grade_mat) <- rownames(RoE_fix)
  colnames(grade_mat) <- colnames(RoE_fix)
  for (i in seq_len(nrow(RoE_fix))) {
    for (j in seq_len(ncol(RoE_fix))) {
      val <- RoE_fix[i, j]
      grade_mat[i, j] <- case_when(
        val <= 1 ~ "-",
        val > 1 & val <= 2 ~ "+",
        val > 2 & val <= 4 ~ "++",
        val > 4 ~ "+++"
      )
    }
  }
  
  all_rows <- rownames(RoE_fix)
  disease_rows <- intersect(target_last, all_rows)
  if (group_type == "condition") {
    pathway_rows <- intersect(fixed_pathway_order, all_rows)
    row_order <- c(pathway_rows, disease_rows)
  } else {
    gene_path_df <- distinct(meta, gene2, pathway) %>% filter(!gene2 %in% target_last)
    gene_path_df$pathway <- factor(gene_path_df$pathway, levels = fixed_pathway_order)
    sorted_gene <- gene_path_df %>% arrange(pathway, gene2) %>% pull(gene2)
    row_order <- c(sorted_gene, disease_rows)
  }
  
  RoE_final <- RoE_fix[row_order, ]
  Freq_raw_final <- Freq_raw_fix[row_order, ]
  Freq_norm_final <- Freq_norm_fix[row_order, ]
  grade_final <- grade_mat[row_order, ]
  
  df_roe <- as.data.frame(RoE_final) %>%
    rownames_to_column(row_var) %>%
    pivot_longer(-all_of(row_var), names_to = "celltype", values_to = "RoE")
  
  df_raw <- as.data.frame(Freq_raw_final) %>%
    rownames_to_column(row_var) %>%
    pivot_longer(-all_of(row_var), names_to = "celltype", values_to = "Freq_raw")
  
  df_norm <- as.data.frame(Freq_norm_final) %>%
    rownames_to_column(row_var) %>%
    pivot_longer(-all_of(row_var), names_to = "celltype", values_to = "Freq_norm")
  
  plot_df <- df_roe %>%
    inner_join(df_raw, by = c(row_var, "celltype")) %>%
    inner_join(df_norm, by = c(row_var, "celltype"))
  
  plot_df <- plot_df %>%
    mutate(
      roe_level = case_when(
        RoE <= 1 ~ "0-1",
        RoE > 1 & RoE <= 2 ~ "1-2",
        RoE > 2 & RoE <= 4 ~ "2-4",
        RoE > 4 ~ ">4"
      ),
      roe_level = factor(roe_level, levels = c("0-1", "1-2", "2-4", ">4")),
      size_grade = case_when(
        roe_level == "0-1"  ~ base_point_size * 0.3,
        roe_level == "1-2"  ~ base_point_size * 0.6,
        roe_level == "2-4"  ~ base_point_size * 0.9,
        roe_level == ">4"   ~ base_point_size * 1.2
      )
    )
  
  plot_df$celltype <- factor(plot_df$celltype, levels = cell_order)
  plot_df[[row_var]] <- factor(plot_df[[row_var]], levels = row_order)
  
  if (freq_mode == "count") {
    fill_col <- "Freq_raw"
    fill_lab <- "Cell proportion(raw)"
  } else {
    fill_col <- "Freq_norm"
    fill_lab <- "Cell proportion(normalized 0-1)"
  }
  
  x_text_theme <- element_text(angle = 45, hjust = 1, colour = "black")
  
  p <- ggplot(plot_df, aes(x = celltype, y = .data[[row_var]])) +
    geom_point(aes(size = size_grade, fill = .data[[fill_col]]), shape = 22, stroke = 0) +
    scale_color_gradientn(colours = prop_color) +
    scale_fill_gradientn(colours = prop_color) +
    scale_size(
      name = "Ro/e level",
      breaks = c(base_point_size*0.3, base_point_size*0.6, base_point_size*0.9, base_point_size*1.2),
      labels = c("0-1", "1-2", "2-4", ">4"),
      range = range(plot_df$size_grade)
    ) +
    labs(fill = fill_lab, x = NULL, y = NULL) +
    guides(
      fill = guide_colorbar(barwidth = legend_key_size, barheight = legend_key_size * 5),
      size = guide_legend(keywidth = legend_key_size, keyheight = legend_key_size, reverse = TRUE)
    ) +
    theme_classic() +
    theme(
      axis.text.x = x_text_theme,
      axis.text.y = element_text(colour = "black"),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      legend.text = element_text(colour = "black"),
      legend.title = element_text(colour = "black")
    )
  
  res_list <- list(
    figure = p,
    RoE_matrix = RoE_final,
    Prop_raw_matrix = Freq_raw_final,
    Prop_norm_matrix = Freq_norm_final
  )
  return(res_list)
}



out_condition_count <- plot_RoE_prop_square_bubble(
  seurat_obj = pbmc1,
  gene_list = gene_list,
  group_type = "condition",
  freq_mode = "count",
  base_point_size = 6,
  legend_key_size = 0.5,
  color_vec =  pal_orrd_100
)
out_condition_count$figure

######2 add 4 forms: original counts, min-max normalized; and subgroup divided

#' Square bubble plot with RoE enrichment and multi-mode cell proportion mapping
#' @param seurat_obj Seurat object, meta must have Level1,Level2,samples,gene columns
#' @param gene_list Named list mapping pathways/genes to gene vectors
#' @param group_type "condition"(disease grouping) or "gene"(pathway grouping)
#' @param freq_mode global_raw/global_norm/subgroup_raw/subgroup_norm; raw uses auto color range, norm fixed 0-1
#' @param base_point_size Base square tile size, scaled by RoE value
#' @param legend_key_size Legend dimension scaling factor
#' @param color_vec Custom fill gradient palette; default blue palette if NULL
#' @return List: figure(ggplot), RoE_matrix, global_raw_matrix, global_norm_matrix, subgroup_raw_matrix, subgroup_norm_matrix, RoE_grade
#' @details Dynamic cell subtype order from Level1 hierarchy; tile fill = proportion, tile size = RoE fold change; X labels truncated to text before "_"

library(scales)
pal_orrd <- brewer.pal(n = 9, name = "OrRd")
pal_orrd_100 <- colorRampPalette(pal_orrd)(1000)
pal_orrd <- rev(brewer.pal(n = 9, name = "RdGy"))
pal_orrd_100 <- colorRampPalette(pal_orrd)(101)
pal_orrd_100<-pal_orrd_100[30:100]

library(scales);library(patchwork)
# pal_orrd <- brewer.pal(n=9,name="OrRd");pal_orrd_100 <- colorRampPalette(pal_orrd)(100)
# rdbu_brewer <- brewer.pal(n=11,name="RdBu");rdbu_100 <- colorRampPalette(rdbu_brewer)(100)
my_colors <- c("Naive CD4 T"="#847AB3","Non-naive CD4 T"="#C8BBD5","Naive CD8 T"="#F6E1EE","Non-naive CD8 T"="#0B71AB","NK"="#96AF95","B cell"="#003D81","Monocyte"="#F49D5C","LDG"="#C95968","Plasma"="#B3D1E7","Erythrocyte"="#FFEFC1","Platelet"="#B3928B","UTC"="#86C7B4","DC"="#F0E442","pDC"="#B07A99","Basophil"="#8C2522")
fixed_level1_order <- c("DC","Monocyte","pDC","LDG","Basophil","B cell","Plasma","Naive CD4 T","Non-naive CD4 T","Naive CD8 T","Non-naive CD8 T","UTC","NK","Erythrocyte","Platelet")
level1_color_map <- my_colors[fixed_level1_order];names(level1_color_map) <- fixed_level1_order

plot_RoE_prop_square_bubble <- function(seurat_obj,gene_list,group_type=c("condition","gene"),
freq_mode=c("global_raw","global_norm","subgroup_raw","subgroup_norm"),base_point_size=5,
legend_key_size=0.5,color_vec=NULL){
  library(dplyr);library(stringr);library(tidyverse);library(RColorBrewer);library(patchwork)
  group_type<-match.arg(group_type);freq_mode<-match.arg(freq_mode)
  meta <- seurat_obj@meta.data
  meta$Level1 <- factor(meta$Level1,levels=fixed_level1_order)
  level2_order <- meta %>% distinct(Level1,Level2) %>% arrange(Level1,Level2) %>% pull(Level2) %>% unique()
  meta$Level2 <- factor(meta$Level2,levels=level2_order)
  meta <- meta %>% arrange(Level1,Level2)
  target_last <- c("SLE","JIA","KD","HC")
  fixed_pathway_order <- c("NFKB pathway","Osteoclast function","Ca2+_flux-PLC","Endolysosomal_nucleic_acid_sensing","Arachidonic acid metabolism","Inflammasome_IL-1β","Immune_metabolic","Inborn_errors_of_cell_death","Negative_regulation_of_IFN-I","Cytoskeleton_and_small_GTPase","TBK1_IRF3","Protein_homeostasis","Uncategorized")
  prop_color <- if(is.null(color_vec)) colorRampPalette(c("#E6F7E6","#B2E8D8","#80D8D8","#3388BB","#0055BB"))(100) else colorRampPalette(color_vec)(100)
  if(group_type=="condition"){
    gene_map <- stack(gene_list);colnames(gene_map) <- c("gene","condition0")
    meta$condition0 <- gene_map$condition0[match(meta$gene,gene_map$gene)]
    meta <- meta %>% mutate(condition=case_when(str_detect(samples,"SLE")~"SLE",str_detect(samples,"JIA")~"JIA",str_detect(samples,"KD")~"KD",str_detect(samples,"HC")~"HC",TRUE~condition0))
    tbl <- table(meta$condition,meta$Level2);row_var <- "condition"
  }else{
    gene_map <- stack(gene_list);colnames(gene_map) <- c("gene","pathway")
    meta$pathway <- gene_map$pathway[match(meta$gene,gene_map$gene)]
    meta <- meta %>% mutate(gene2=gene) %>% mutate(gene2=case_when(str_detect(samples,"SLE")~"SLE",str_detect(samples,"JIA")~"JIA",str_detect(samples,"KD")~"KD",str_detect(samples,"HC")~"HC",TRUE~gene2))
    tbl <- table(meta$gene2,meta$Level2);row_var <- "gene2"
  }
  cell_custom_order <- levels(meta$Level2);cell_order <- intersect(cell_custom_order,colnames(tbl))
  obs_mat <- as.data.frame.matrix(tbl);obs_mat <- obs_mat[,match(cell_order,colnames(obs_mat))]
  total <- sum(obs_mat);row_sum <- rowSums(obs_mat);col_sum <- colSums(obs_mat)
  E <- outer(row_sum,col_sum)/total;RoE <- obs_mat/E;RoE_fix <- RoE[,cell_order]
  freq <- as.data.frame(t(apply(obs_mat,1,function(x)as.numeric(x/sum(x)))));colnames(freq) <- colnames(obs_mat)
  freq_norm <- as.data.frame(apply(freq,2,function(x)(x-min(x))/(max(x)-min(x))))
  allsubcol <- grep("_",colnames(obs_mat),value=T);allsubcol <- sapply(allsubcol,function(x)unlist(strsplit(x,"_"))[1])
  subcol <- c("DC0","Mono0","B0","Plasma0","NnCD4T0","NnCD8T0","UTC0","NK0");subgroup_freq <- data.frame()
  for(i in 1:nrow(obs_mat)){
    y <- obs_mat[i,];y1 <- data.frame(value=as.numeric(y),cells=colnames(y))
    y1$subcol <- sapply(y1$cells,function(x){hit<-subcol[sapply(subcol,function(z)grepl(paste0("^",z),x))];if(length(hit)>0)hit[1]else"total"})
    y1 <- y1 %>% group_by(subcol) %>% mutate(total_groups=sum(value),prop=value/total_groups)
    y2 <- as.numeric(y1$prop);names(y2) <- y1$cells;subgroup_freq <- rbind(subgroup_freq,y2)
  }
  colnames(subgroup_freq) <- colnames(obs_mat);rownames(subgroup_freq) <- rownames(obs_mat)
  subgroup_freq_norm <- as.data.frame(apply(subgroup_freq,2,function(x)(x-min(x))/(max(x)-min(x))))
  grade_mat <- matrix("",nrow=nrow(RoE_fix),ncol=ncol(RoE_fix),dimnames=dimnames(RoE_fix))
  for(i in seq_len(nrow(RoE_fix))){
    for(j in seq_len(ncol(RoE_fix))){
      v <- RoE_fix[i,j];grade_mat[i,j] <- case_when(v<=1~"-",v>1&v<=2~"+",v>2&v<=4~"++",TRUE~"+++")
    }
  }
  all_rows <- rownames(RoE_fix);disease_rows <- intersect(target_last,all_rows)
  if(group_type=="condition"){
    pathway_rows <- intersect(fixed_pathway_order,all_rows);row_order <- c(pathway_rows,disease_rows)
  }else{
    gene_path_df <- distinct(meta,gene2,pathway) %>% filter(!gene2%in%target_last)
    gene_path_df$pathway <- factor(gene_path_df$pathway,levels=fixed_pathway_order)
    sorted_gene <- gene_path_df %>% arrange(pathway,gene2) %>% pull(gene2)
    row_order <- c(sorted_gene,disease_rows)
  }
  RoE_final <- RoE_fix[row_order,];freq_final <- freq[row_order,cell_order]
  freq_norm_final <- freq_norm[row_order,cell_order];subgroup_freq_final <- subgroup_freq[row_order,cell_order]
  subgroup_freq_norm_final <- subgroup_freq_norm[row_order,cell_order];grade_final <- grade_mat[row_order,cell_order]
  plot_mat <- switch(freq_mode,global_raw=freq_final,global_norm=freq_norm_final,subgroup_raw=subgroup_freq_final,subgroup_norm=subgroup_freq_norm_final)
  fill_lab <- switch(freq_mode,global_raw="Raw subgroup fraction",global_norm="Normalized fraction (0-1)",subgroup_raw="Raw cell fraction per group",subgroup_norm="Subgroup normalized fraction (0-1)")
  short_cell_label <- function(nm)if(str_detect(nm,"_"))str_split(nm,"_")[[1]][1]else nm
  df_prop <- as.data.frame(plot_mat) %>% rownames_to_column(row_var) %>% pivot_longer(-all_of(row_var),names_to="cell_full",values_to="prop_val") %>% mutate(cell_short=sapply(cell_full,short_cell_label))
  df_roe <- as.data.frame(RoE_final) %>% rownames_to_column(row_var) %>% pivot_longer(-all_of(row_var),names_to="cell_full",values_to="RoE")
  df_grade <- as.data.frame(grade_final) %>% rownames_to_column(row_var) %>% pivot_longer(-all_of(row_var),names_to="cell_full",values_to="roe_grade")
  plot_df <- df_prop %>% inner_join(df_roe,by=c(row_var,"cell_full")) %>% inner_join(df_grade,by=c(row_var,"cell_full"))
  plot_df <- plot_df %>% mutate(roe_level=case_when(RoE<=1~"0-1",RoE>1&RoE<=2~"1-2",RoE>2&RoE<=4~"2-4",RoE>4~">4"),roe_level=factor(roe_level,levels=c("0-1","1-2","2-4",">4")),size_grade=case_when(roe_level=="0-1"~base_point_size*0.3*0.75,roe_level=="1-2"~base_point_size*0.6*0.75,roe_level=="2-4"~base_point_size*0.9*0.75,roe_level==">4"~base_point_size*1.2*0.75))
  # Low RoE stat table (RoE <= 1)
  low_roe_df <- plot_df %>% filter(RoE <= 1)
  cell_level1_map <- meta %>% distinct(Level2,Level1)
  low_roe_df <- low_roe_df %>% left_join(cell_level1_map,by=c("cell_full"="Level2"))
  stat_low_roe <- low_roe_df %>% group_by(.data[[row_var]]) %>% summarise(count_lowRoE=n(),cell_lineage=paste(unique(Level1),collapse=", "),.groups="drop")
  colnames(stat_low_roe) <- c("group","count_lowRoE","cell_lineage")
  # Filter plot data, remove RoE<=1 subtype
  plot_df_filter <- plot_df %>% filter(roe_level != "0-1")
  plot_df_filter$cell_full <- factor(plot_df_filter$cell_full,levels=cell_order)
  plot_df_filter[[row_var]] <- factor(plot_df_filter[[row_var]],levels=row_order)
  plot_df_filter <- plot_df_filter %>% left_join(cell_level1_map,by=c("cell_full"="Level2"))
  plot_df_filter$Level1 <- factor(plot_df_filter$Level1,levels=fixed_level1_order)
  x_text_theme <- element_text(angle=90,hjust=1,vjust=0.3,colour="black")
  p_main <- ggplot(plot_df_filter,aes(x=cell_full,y=.data[[row_var]])) + geom_point(aes(size=size_grade,fill=prop_val),shape=22,stroke=0) + scale_fill_gradientn(colours=prop_color,limits=range(plot_df_filter$prop_val,na.rm=T)) + scale_size(name="RoE level",breaks=c(base_point_size*0.3*0.75,base_point_size*0.6*0.75,base_point_size*0.9*0.75,base_point_size*1.2*0.75),labels=c("0-1","1-2","2-4",">4"),range=range(plot_df_filter$size_grade)) + labs(fill=fill_lab,x=NULL,y=NULL) + guides(fill=guide_colorbar(barwidth=legend_key_size*8,barheight=legend_key_size),size=guide_legend(keywidth=legend_key_size,keyheight=legend_key_size,reverse=T)) + theme_classic() + theme(axis.text.x=x_text_theme,axis.text.y=element_text(colour="black"),axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),panel.border=element_blank(),axis.line=element_blank(),legend.text=element_text(colour="black"),legend.title=element_text(colour="black"),legend.position="top",legend.box="horizontal") + scale_x_discrete(labels=function(x)sapply(x,short_cell_label))
  # Prepare percentage fill bar data, fix undefined stat_full
  stat_raw <- plot_df_filter %>% filter(RoE > 2) %>% group_by(.data[[row_var]],Level1) %>% summarise(count=n(),.groups="drop")
  all_row_names <- unique(stat_raw[[row_var]])
  full_grid <- expand.grid(temp_row=all_row_names,Level1=fixed_level1_order)
  colnames(full_grid)[1] <- row_var
  stat_full <- full_grid %>% left_join(stat_raw,by=c(row_var,"Level1")) %>% mutate(count=replace_na(count,0))
  stat_df <- stat_full %>% group_by(.data[[row_var]]) %>% mutate(percent=count/sum(count)) %>% ungroup()
  stat_df[[row_var]] <- factor(stat_df[[row_var]],levels=row_order)
  stat_df$Level1 <- factor(stat_df$Level1,levels=rev(fixed_level1_order))
  p_bar <- ggplot(stat_df,aes(y=.data[[row_var]],x=percent,fill=Level1)) + geom_col(position="fill",width=0.7) + scale_fill_manual(values=rev(level1_color_map),limits=rev(fixed_level1_order),labels=fixed_level1_order) + scale_x_continuous(labels=percent_format(),expand=c(0,0)) + labs(x="Proportion of subtypes with RoE > 2",y=NULL) + theme_classic() + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position="none",axis.line.y=element_blank())
  p_combine <- p_main + p_bar + plot_layout(widths=c(9,1))
  return(list(full_plot=p_combine,
              bubble_plot=p_main,
              percent_bar=p_bar,
              lowRoE_summary=stat_low_roe,
              raw_freq=freq_final,
              raw_freqnorm=freq_norm_final,
              subgroup_freq=subgroup_freq_final,
              subgroup_freqnorm=subgroup_freq_norm_final))
}

res <- plot_RoE_prop_square_bubble(seurat_obj=pbmc1,gene_list=gene_list,
group_type="condition",freq_mode="subgroup_norm",base_point_size=5,
legend_key_size=0.47,color_vec=pal_orrd_100)

res$full_plot
res$lowRoE_summary

#####3 cell prop for each sample at Level2, and compare with different status (HC MP PP)
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggsignif)

status_order <- c("unaffected control", "polygenic patient", "monogenic patient")
status_col <- c("unaffected control"="#4375A0","polygenic patient"="#B18749","monogenic patient"="#912623")
subcol <- c("DC0","Mono0","B0","Plasma0","NnCD4T0","NnCD8T0","UTC0","NK0")
fixed_level1_order <- c("DC","Monocyte","pDC","LDG","Basophil","B cell","Plasma","Naive CD4 T","Non-naive CD4 T","Naive CD8 T","Non-naive CD8 T","UTC","NK","Erythrocyte","Platelet")
dodge_pos <- position_dodge(width = 0.6)

meta <- pbmc1@meta.data
meta$Level1 <- factor(meta$Level1, levels = fixed_level1_order)
level2_order <- meta %>% distinct(Level1, Level2) %>% arrange(Level1, Level2) %>% pull(Level2) %>% unique()
meta$Level2 <- factor(meta$Level2, levels = level2_order)
meta <- meta %>% arrange(Level1, Level2)

tbl <- table(meta$samples, meta$Level2)
obs_mat <- as.data.frame.matrix(tbl)
cell_custom_order <- levels(meta$Level2)
cell_order <- intersect(cell_custom_order, colnames(tbl))
obs_mat <- obs_mat[, match(cell_order, colnames(obs_mat))]

subgroup_freq <- data.frame()
for (i in 1:nrow(obs_mat)) {
  y <- obs_mat[i, ]
  y1 <- data.frame(value = as.numeric(y), cells = colnames(y))
  y1$subcol <- sapply(y1$cells, function(x) {
    hit <- subcol[sapply(subcol, function(z) grepl(paste0("^", z), x))]
    if (length(hit) > 0) hit[1] else "total"
  })
  y1 <- y1 %>% group_by(subcol) %>% mutate(total_groups = sum(value))
  y1 <- y1 %>% mutate(prop = case_when(total_groups == 0 ~ 0, TRUE ~ value / total_groups))
  y2 <- as.numeric(y1$prop)
  names(y2) <- y1$cells
  subgroup_freq <- rbind(subgroup_freq, y2)
}
sum(is.na(subgroup_freq));min(subgroup_freq)

#subgroup_freq[is.na(subgroup_freq)] <- 0
colnames(subgroup_freq) <- colnames(obs_mat)
rownames(subgroup_freq) <- rownames(obs_mat)

subgroup_freq_norm <- apply(subgroup_freq, 2, function(x) {
  diff_val <- max(x) - min(x)
  if (diff_val == 0) return(rep(0, length(x)))
  res <- (x - min(x)) / diff_val
  res[is.na(res) | is.infinite(res)] <- 0
  res
}) %>% as.data.frame()
rownames(subgroup_freq_norm) <- rownames(subgroup_freq_norm)
colnames(subgroup_freq_norm) <- colnames(subgroup_freq_norm)
sum(is.na(subgroup_freq_norm));min(subgroup_freq_norm)


raw_long <- subgroup_freq %>% rownames_to_column("samples") %>% pivot_longer(-samples, names_to = "Level2", values_to = "raw_prop")
norm_long <- subgroup_freq_norm %>% rownames_to_column("samples") %>% pivot_longer(-samples, names_to = "Level2", values_to = "norm_prop")
sample_info <- meta %>% distinct(samples, status)
plot_df <- raw_long %>%
  left_join(norm_long, by = c("samples", "Level2")) %>%
  left_join(sample_info, by = "samples") %>%
  mutate(status = factor(status, levels = status_order), Level2 = factor(Level2, levels = level2_order))
sum(is.na(plot_df$raw_prop))
sum(is.na(plot_df$norm_prop))
dim(plot_df)
plot_df <- plot_df %>% filter(!is.na(raw_prop), !is.infinite(raw_prop), !is.na(norm_prop), !is.infinite(norm_prop))
dim(plot_df)

kruskal_safe <- function(df, val_col) {
  g <- unique(df$status)
  if (length(g) < 2) return(tibble(p_value = NA_real_))
  form <- as.formula(paste0(val_col, " ~ status"))
  tibble(p_value = kruskal.test(form, data = df)$p.value)
}

level2_group_info <- plot_df %>% group_by(Level2) %>% summarise(group_total = n_distinct(status), .groups = "drop")
table(level2_group_info$group_total)

# 5. Kruskal test + Bonferroni correction: raw proportion
kw_raw_all <- plot_df %>%
  group_by(Level2) %>%
  group_modify(~kruskal_safe(.x, "raw_prop")) %>%
  ungroup() %>%
  left_join(level2_group_info, by = "Level2")
sum(is.na(kw_raw_all$p_value))

valid_p_raw <- !is.na(kw_raw_all$p_value)
kw_raw_all$padj <- NA_real_
kw_raw_all$padj[valid_p_raw] <- p.adjust(kw_raw_all$p_value[valid_p_raw], method = "bonferroni")
sum(is.na(kw_raw_all$padj))

kw_raw_all <- kw_raw_all %>% mutate(
  signif = case_when(
    group_total < 2 ~ "single_group",
    padj < 0.001 ~ "***",
    padj < 0.01 ~ "**",
    padj < 0.05 ~ "*",
    TRUE ~ "ns"
  ),
  sig_flag = ifelse(group_total >= 2 & padj < 0.05, "sig", "non_sig")
) %>% as.data.frame()

kw_norm_all <- plot_df %>%
  group_by(Level2) %>%
  group_modify(~kruskal_safe(.x, "norm_prop")) %>%
  ungroup() %>%
  left_join(level2_group_info, by = "Level2")

sum(is.na(kw_norm_all$p_value))
valid_p_norm <- !is.na(kw_norm_all$p_value)
kw_norm_all$padj <- NA_real_
kw_norm_all$padj[valid_p_norm] <- p.adjust(kw_norm_all$p_value[valid_p_norm], method = "bonferroni")
sum(is.na(kw_norm_all$padj))

kw_norm_all <- kw_norm_all %>% mutate(
  signif = case_when(
    group_total < 2 ~ "single_group",
    padj < 0.001 ~ "***",
    padj < 0.01 ~ "**",
    padj < 0.05 ~ "*",
    TRUE ~ "ns"
  ),
  sig_flag = ifelse(group_total >= 2 & padj < 0.05, "sig", "non_sig")
) %>% as.data.frame()

sig_level2_raw <- kw_raw_all %>% filter(sig_flag == "sig") %>% pull(Level2)
sig_level2_norm <- kw_norm_all %>% filter(sig_flag == "sig") %>% pull(Level2)

plot_raw_sig <- plot_df %>%
  filter(Level2 %in% sig_level2_raw) %>%
  mutate(Level2 = factor(Level2, levels = intersect(level2_order, sig_level2_raw)))

plot_norm_sig <- plot_df %>%
  filter(Level2 %in% sig_level2_norm) %>%
  mutate(Level2 = factor(Level2, levels = intersect(level2_order, sig_level2_norm)))

short_label_fun <- function(x) {
  str_replace_all(x, "^(DC|B|Mono|Plasma|NnCD4T|NnCD8T|UTC|NK)(\\d+).*", "\\1\\2")
}

#(1) Boxplot: min-max normalized proportion, only significant features
box_raw <- ggplot(plot_raw_sig, aes(x = Level2, y = raw_prop, color = status)) +
  geom_boxplot(position = dodge_pos, fill = NA, width = 0.55, linewidth = 0.6, outlier.shape = NA) +
  scale_color_manual(values = status_col) +
  labs(x = "", y = "Raw subgroup proportion (within parent Level1)") +
  scale_x_discrete(labels = short_label_fun) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 35, hjust = 1, colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line("black"),
    panel.border = element_blank(),
    axis.line = element_line("black"))

#(2) Boxplot: min-max normalized proportion, only significant features and remove too lower values
library(dplyr)
library(ggplot2)
library(stringr)

# label function
short_label_fun <- function(x) {
  str_extract(x, "^[^_]+")
}

box_norm <- ggplot(
  plot_norm_sig[!plot_norm_sig$Level2 %in% c("DC05_cDC2_LY96","NnCD4T05_EMRA_GNLY"),], 
  aes(x = Level2, y = norm_prop, color = status)
) +
  geom_boxplot(
    position = position_dodge(width = 0.8), 
    fill = NA, 
    width = 0.7, 
    linewidth = 0.6, 
    outlier.shape = NA
  ) +
  scale_color_manual(values = status_col) +
  labs(x = "", y = "Cell proportion (%,0-1 scale)") +
  scale_x_discrete(labels = short_label_fun) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 35, hjust = 1, colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line("black"),
    panel.border = element_blank(),
    axis.line = element_line("black")
  )

box_norm

#(3) dodge gitter
library(dplyr)
library(ggplot2)
library(stringr)

short_label_fun <- function(x) {
  str_extract(x, "^[^_]+")
}

box_norm <- ggplot(
  plot_norm_sig[!plot_norm_sig$Level2 %in% c("DC05_cDC2_LY96","NnCD4T05_EMRA_GNLY"),], 
  aes(x = Level2, y = norm_prop, color = status)
) +
  geom_boxplot(
    position = position_dodge(width = 0.8), 
    fill = NA, 
    width = 0.7, 
    linewidth = 0.6, 
    outlier.shape = NA # 不绘制箱线自带离群点
  ) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    size = 0.8, # 缩小散点尺寸
    alpha = 0.7
  ) +
  scale_color_manual(values = status_col) +
  labs(x = "", y = "Cell proportion (%,0-1 scale)") +
  scale_x_discrete(labels = short_label_fun) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 35, hjust = 1, colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line("black"),
    panel.border = element_blank(),
    axis.line = element_line("black")
  )

box_norm

#(4) without outier
library(dplyr)
library(ggplot2)
library(stringr)

short_label_fun <- function(x) {
  str_extract(x, "^[^_]+")
}

# Pre-filter data to remove extreme outliers beyond whiskers
plot_clean <- plot_norm_sig %>%
  filter(!Level2 %in% c("DC05_cDC2_LY96","NnCD4T05_EMRA_GNLY")) %>%
  group_by(Level2, status) %>%
  mutate(
    q1 = quantile(norm_prop, 0.25, na.rm = TRUE),
    q3 = quantile(norm_prop, 0.75, na.rm = TRUE),
    iqr = q3 - q1,
    lower_bound = q1 - 1.5 * iqr,
    upper_bound = q3 + 1.5 * iqr
  ) %>%
  filter(norm_prop >= lower_bound & norm_prop <= upper_bound) %>%
  ungroup()

box_norm <- ggplot(
  plot_clean, 
  aes(x = Level2, y = norm_prop, color = status)
) +
  geom_boxplot(
    position = position_dodge(width = 0.8), 
    fill = NA, 
    width = 0.7, 
    linewidth = 0.6, 
    outlier.shape = NA
  ) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    size = 0.8,
    alpha = 0.7
  ) +
  scale_color_manual(values = status_col) +
  labs(x = "", y = "Cell proportion (%,0-1 scale)") +
  scale_x_discrete(
    labels = short_label_fun
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 35, hjust = 1, colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line("black"),
    panel.border = element_blank(),
    axis.line = element_line("black")
  )

box_norm


#(5) facet
library(dplyr)
library(ggplot2)
library(stringr)

plot_norm_filtered <- plot_norm_sig %>% 
  filter(!Level2 %in% c("DC05_cDC2_LY96","NnCD4T05_EMRA_GNLY")) %>%
  mutate(Level2 = str_extract(Level2, "^[^_]+"))

box_norm <- ggplot(plot_norm_filtered, aes(x = status, y = norm_prop, color = status)) +
  geom_boxplot(position = dodge_pos, fill = NA, width = 0.55, linewidth = 0.6,outlier.shape = NA) + # 
  scale_color_manual(values = status_col) +
  labs(x = "", y = "") +
  facet_wrap(~Level2, scales = "free_y", nrow = 3, ncol = 6, drop = TRUE) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
    axis.line = element_blank(),
    strip.background = element_rect(fill = "gray85", colour = "black"),
    strip.text = element_text(colour = "black", face = "bold", size = 8,
                              margin = margin(t = 2, b = 2, l = 4, r = 4)),
    panel.spacing = unit(4, "mm")
  )
box_norm

#(6) facet with gitter
library(dplyr)
library(ggplot2)
library(stringr)

dodge_pos <- position_dodge(width = 0.55)
jit_dodge <- position_jitterdodge(
  dodge.width = 0.55,
  jitter.width = 0.22
)

plot_norm_filtered <- plot_norm_sig %>% 
  filter(!Level2 %in% c("DC05_cDC2_LY96","NnCD4T05_EMRA_GNLY")) %>%
  mutate(Level2 = str_extract(Level2, "^[^_]+"))

box_norm <- ggplot(plot_norm_filtered, aes(x = status, y = norm_prop, color = status)) +
  geom_jitter(position = jit_dodge, size = 0.6) +
  geom_boxplot(
    position = dodge_pos, 
    fill = NA, 
    width = 0.55, 
    linewidth = 0.6,
    outlier.shape = NA
  ) +
  scale_color_manual(values = status_col) +
  labs(x = "", y = "") +
  facet_wrap(~Level2, scales = "free_y", nrow = 3, ncol = 6, drop = TRUE) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
    axis.line = element_blank(),
    strip.background = element_rect(fill = "gray85", colour = "black"),
    strip.text = element_text(colour = "black", face = "bold", size = 8,
                              margin = margin(t = 2, b = 2, l = 4, r = 4)),
    panel.spacing = unit(4, "mm")
  )

box_norm


#####4 cell prop at different age groups
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggsignif)

status_order <- c("unaffected control", "polygenic patient", "monogenic patient")
status_col <- c("unaffected control"="#4375A0","polygenic patient"="#B18749","monogenic patient"="#912623")
subcol <- c("DC0","Mono0","B0","Plasma0","NnCD4T0","NnCD8T0","UTC0","NK0")
fixed_level1_order <- c("DC","Monocyte","pDC","LDG","Basophil","B cell","Plasma","Naive CD4 T","Non-naive CD4 T","Naive CD8 T","Non-naive CD8 T","UTC","NK","Erythrocyte","Platelet")
dodge_pos <- position_dodge(width = 0.6)

meta <- pbmc1@meta.data
meta$Level1 <- factor(meta$Level1, levels = fixed_level1_order)
level2_order <- meta %>% distinct(Level1, Level2) %>% arrange(Level1, Level2) %>% pull(Level2) %>% unique()
meta$Level2 <- factor(meta$Level2, levels = level2_order)
meta <- meta %>% arrange(Level1, Level2)

tbl <- table(meta$samples, meta$Level2)
obs_mat <- as.data.frame.matrix(tbl)
cell_custom_order <- levels(meta$Level2)
cell_order <- intersect(cell_custom_order, colnames(tbl))
obs_mat <- obs_mat[, match(cell_order, colnames(obs_mat))]

subgroup_freq <- data.frame()
for (i in 1:nrow(obs_mat)) {
  y <- obs_mat[i, ]
  y1 <- data.frame(value = as.numeric(y), cells = colnames(y))
  y1$subcol <- sapply(y1$cells, function(x) {
    hit <- subcol[sapply(subcol, function(z) grepl(paste0("^", z), x))]
    if (length(hit) > 0) hit[1] else "total"
  })
  y1 <- y1 %>% group_by(subcol) %>% mutate(total_groups = sum(value))
  y1 <- y1 %>% mutate(prop = case_when(total_groups == 0 ~ 0, TRUE ~ value / total_groups))
  y2 <- as.numeric(y1$prop)
  names(y2) <- y1$cells
  subgroup_freq <- rbind(subgroup_freq, y2)
}
sum(is.na(subgroup_freq));min(subgroup_freq)

#subgroup_freq[is.na(subgroup_freq)] <- 0
colnames(subgroup_freq) <- colnames(obs_mat)
rownames(subgroup_freq) <- rownames(obs_mat)

subgroup_freq_norm <- apply(subgroup_freq, 2, function(x) {
  diff_val <- max(x) - min(x)
  if (diff_val == 0) return(rep(0, length(x)))
  res <- (x - min(x)) / diff_val
  res[is.na(res) | is.infinite(res)] <- 0
  res
}) %>% as.data.frame()
rownames(subgroup_freq_norm) <- rownames(subgroup_freq_norm)
colnames(subgroup_freq_norm) <- colnames(subgroup_freq_norm)
sum(is.na(subgroup_freq_norm));min(subgroup_freq_norm)
dim(subgroup_freq_norm)

metadata<-as.data.frame(read.xlsx("E:\\AID cohort\\code\\result\\metadata.xlsx"))
metadata$source <- ifelse(
  grepl("^in-", metadata$resource),
  "inhouse",
  "public"
)
unique(metadata$source)
unique(metadata$age)
idx_m <- grepl("m", metadata$age, ignore.case = TRUE)
metadata$age_year[idx_m] <- as.numeric(gsub("[^0-9.]", "", metadata$age[idx_m])) / 12
idx_y <- grepl("y", metadata$age, ignore.case = TRUE)
metadata$age_year[idx_y] <- as.numeric(gsub("[^0-9.]", "", metadata$age[idx_y]))
metadata$age_year
m2<-metadata[-which(is.na(metadata$age_year)),c("dataset","age_year","status")]

subgroup_freq[match(m2$dataset,rownames(subgroup_freq_norm)),]->agedata
agedata$samples<-rownames(agedata)
#identical(agedata$samples,m2$dataset)
agedata$age<-m2$age_year
agedata$status<-m2$status

## (1) Without healthy control calculation
agedata <- agedata %>%
  mutate(
    age_group = case_when(
      age < 4 ~ "<4",
      age >= 4 & age < 18 ~ "4–18",
      age >= 18 ~ "≥18"
    )
  )

sample_count <- agedata %>%
  count(status, age_group, name = "n_sample") %>%
  pivot_wider(names_from = age_group, values_from = n_sample, values_fill = 0)
print(sample_count)

df_all <- agedata %>%
  filter(status %in% c("monogenic patient", "polygenic patient")) %>%
  mutate(
    age_group = case_when(
      age < 4 ~ "<4",
      age >= 4 & age < 18 ~ "4–18",
      age >= 18 ~ "≥18"
    ),
    age_group = factor(age_group, levels = c("<4", "4–18", "≥18"))
  )

cell_cols <- setdiff(colnames(df_all), c("samples", "age", "status", "age_group"))

long_df <- df_all %>%
  pivot_longer(
    cols = all_of(cell_cols),
    names_to = "cell",
    values_to = "prop"
  )

get_age_kruskal_res <- function(target_status) {
  sub <- long_df %>% filter(status == target_status)
  
  mean_tab <- sub %>%
    group_by(cell, age_group) %>%
    summarise(mean_prop = mean(prop), .groups = "drop") %>%
    pivot_wider(names_from = age_group, values_from = mean_prop)
  
  p_tab <- sub %>%
    group_by(cell) %>%
    group_modify(~{
      kt <- kruskal.test(prop ~ age_group, data = .x)
      tibble(p_raw = kt$p.value)
    }) %>% ungroup()
  
  res_full <- inner_join(mean_tab, p_tab, by = "cell") %>%
    mutate(p_adj = p.adjust(p_raw, method = "BH")) %>%
    select(cell, `<4`, `4–18`, `≥18`, p_raw, p_adj)
  
  return(res_full)
}

res_mono <- as.data.frame(get_age_kruskal_res("monogenic patient"))
res_poly <- as.data.frame(get_age_kruskal_res("polygenic patient"))

head(res_mono)
res_mono[res_mono$p_adj < 0.05, ]
head(res_poly)
res_poly[res_poly$p_adj < 0.05, ]


## (2) With healthy control evaluation & log2FC calculation
agedata <- agedata %>%
  mutate(
    age_group = case_when(
      age < 4 ~ "lt4",
      age >= 4 & age < 18 ~ "mid",
      age >= 18 ~ "ge18"
    ),
    age_group = factor(age_group, levels = c("lt4", "mid", "ge18"))
  )

sample_count <- agedata %>%
  count(status, age_group, name = "n_sample") %>%
  pivot_wider(names_from = age_group, values_from = n_sample, values_fill = 0)
print(sample_count)

cell_cols <- setdiff(colnames(agedata), c("samples", "age", "status", "age_group"))

long_all <- agedata %>%
  pivot_longer(
    cols = all_of(cell_cols),
    names_to = "cell",
    values_to = "prop"
  )

all_age_lvls <- c("lt4", "mid", "ge18")
all_cells <- unique(long_all$cell)

hc_raw <- long_all %>%
  filter(status == "unaffected control") %>%
  group_by(cell, age_group) %>%
  summarise(mean_val = mean(prop), .groups = "drop")

hc_grid <- expand.grid(cell = all_cells, age_group = all_age_lvls)
hc_tab <- hc_grid %>%
  left_join(hc_raw, by = c("cell", "age_group")) %>%
  mutate(mean_val = replace_na(mean_val, 0)) %>%
  pivot_wider(
    names_from = age_group,
    values_from = mean_val,
    names_glue = "ctrl_{age_group}"
  )

get_sig_result <- function(target_status) {
  sub_pat <- long_all %>% filter(status == target_status)
  
  pat_raw <- sub_pat %>%
    group_by(cell, age_group) %>%
    summarise(mean_val = mean(prop), .groups = "drop")
  
  pat_grid <- expand.grid(cell = unique(sub_pat$cell), age_group = all_age_lvls)
  pat_tab <- pat_grid %>%
    left_join(pat_raw, by = c("cell", "age_group")) %>%
    mutate(mean_val = replace_na(mean_val, 0)) %>%
    pivot_wider(
      names_from = age_group,
      values_from = mean_val,
      names_glue = "pat_{age_group}"
    )
  
  p_tab <- sub_pat %>%
    group_by(cell) %>%
    group_modify(~tibble(p_raw = kruskal.test(prop ~ age_group, data = .x)$p.value)) %>%
    ungroup()
  
  res_full <- pat_tab %>%
    inner_join(hc_tab, by = "cell") %>%
    inner_join(p_tab, by = "cell") %>%
    mutate(
      p_adj = p.adjust(p_raw, method = "BH"),
      log2FC_lt4  = log2((pat_lt4  + 0.01) / (ctrl_lt4  + 0.01)),
      log2FC_mid  = log2((pat_mid  + 0.01) / (ctrl_mid  + 0.01)),
      log2FC_ge18 = log2((pat_ge18 + 0.01) / (ctrl_ge18 + 0.01))
    ) %>%
    filter(!is.na(p_adj) & !is.nan(p_adj) & p_adj < 0.05) %>%
    select(
      cell,
      pat_lt4, pat_mid, pat_ge18,
      ctrl_lt4, ctrl_mid, ctrl_ge18,
      log2FC_lt4, log2FC_mid, log2FC_ge18,
      p_raw, p_adj
    )
  
  return(as.data.frame(res_full))
}

res_mono_sig <- get_sig_result("monogenic patient")
res_poly_sig <- get_sig_result("polygenic patient")

head(res_mono_sig)
head(res_poly_sig)

##(3) bubble plot
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

res_mono_sig <- res_mono_sig %>% mutate(group = "monogenic")
res_poly_sig <- res_poly_sig %>% mutate(group = "polygenic")
plot_df <- bind_rows(res_mono_sig, res_poly_sig)

plot_df <- plot_df %>% mutate(cell_short = case_when(str_detect(cell, "_") ~ str_split(cell, "_", simplify = TRUE)[,1], TRUE ~ cell))
mono_short <- unique(str_split(res_mono_sig$cell, "_", simplify = TRUE)[,1])
poly_short <- unique(str_split(res_poly_sig$cell, "_", simplify = TRUE)[,1])
cell_order <- c(mono_short, poly_short)

plot_long <- plot_df %>% pivot_longer(cols = c(pat_lt4, pat_mid, pat_ge18, log2FC_lt4, log2FC_mid, log2FC_ge18), names_to = c(".value", "age"), names_sep = "_") %>% mutate(age = factor(age, levels = c("lt4", "mid", "ge18"), labels = c("<4", "4–18", "≥18")), cell_short = factor(cell_short, levels = cell_order))
bar_df <- plot_df %>% distinct(cell_short, p_adj) %>% mutate(cell_short = factor(cell_short, levels = cell_order), neg_log10_padj = -log10(p_adj), threshold = -log10(0.05))

custom_pal <- colorRampPalette(c("#1c2b5a","#192a57","#4E70AF","#7091C7","#9EBCDB","#C8D6E7","#E8EDF1","#F2EBE5","#ECD0B4","#C16D58","#A13D3B","#951d22"))(101)

p_bar <- ggplot(bar_df, aes(x = cell_short, y = neg_log10_padj)) +
  geom_col(fill = "#847AB3", width = 0.6) +
  geom_hline(yintercept = unique(bar_df$threshold), linetype = "dashed", colour = "#888888", linewidth = 0.3) +
  labs(x = NULL, y = "") +
  theme_classic() +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.text.y = element_text(colour = "black"),axis.title.y = element_text(colour = "black"),axis.ticks.y = element_line(colour = "black"),axis.line.y = element_line(colour = "black"),panel.grid = element_blank(),panel.border = element_blank())

p_bubble <- ggplot(plot_long, aes(x = cell_short, y = age)) +
  geom_point(aes(size = pat, fill = log2FC), shape = 21, stroke = 0.3) +
  scale_fill_gradientn(colours = custom_pal, name = NULL) +
  scale_size(range = c(2, 12), name = NULL) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, 
  colour = "black"),axis.text.y = element_text(colour = "black"),
  axis.ticks = element_line(colour = "black"),
  panel.grid.major = element_line(colour = "#cccccc", linewidth = 0.3,
  linetype = "dashed"),panel.grid.minor = element_line(colour = "#e6e6e6", linewidth = 0.2, linetype = "dashed"),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),axis.line = element_blank())
p_comb <- p_bar / p_bubble + plot_layout(heights = c(1, 3))
print(p_comb)
