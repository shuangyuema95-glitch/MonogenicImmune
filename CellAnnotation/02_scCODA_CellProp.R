library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(tidyverse)
library(stringr)
library(patchwork)
library(cowplot)

setwd("E:\\AID cohort\\code\\result")
pbmc1<-readRDS("PBMC1.rds")

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

fixed_level1_order <- c("DC","Monocyte","pDC","LDG","Basophil","B cell","Plasma","Naive CD4 T","Non-naive CD4 T","Naive CD8 T","Non-naive CD8 T","UTC","NK","Erythrocyte","Platelet")


######1 count matrix input
make_scCODA_count_matrix <- function(
    seurat_obj,
    group_type = c("condition", "gene"),
    gene_list,
    fixed_level1_order
){
  group_type <- match.arg(group_type)
  meta <- seurat_obj@meta.data
  meta$Level1 <- factor(meta$Level1,levels=fixed_level1_order)
  level2_order <- meta %>% distinct(Level1,Level2) %>% arrange(Level1,Level2) %>% pull(Level2) %>% unique()
  meta$Level2 <- factor(meta$Level2,levels=level2_order)
  meta <- meta %>% arrange(Level1,Level2)
  
  if(group_type=="condition"){
    gene_map <- stack(gene_list)
    colnames(gene_map) <- c("gene","condition0")
    meta$condition0 <- gene_map$condition0[match(meta$gene,gene_map$gene)]
    meta <- meta %>% mutate(condition=case_when(
      str_detect(samples,"SLE")~"SLE",
      str_detect(samples,"JIA")~"JIA",
      str_detect(samples,"KD")~"KD",
      str_detect(samples,"HC")~"HC",
      TRUE~condition0
    ))
    tbl <- table(meta$condition,meta$Level2)
    row_var <- "condition"
  }else{
    gene_map <- stack(gene_list)
    colnames(gene_map) <- c("gene","pathway")
    meta$pathway <- gene_map$pathway[match(meta$gene,gene_map$gene)]
    meta <- meta %>% mutate(gene2=gene) %>% mutate(gene2=case_when(
      str_detect(samples,"SLE")~"SLE",
      str_detect(samples,"JIA")~"JIA",
      str_detect(samples,"KD")~"KD",
      str_detect(samples,"HC")~"HC",
      TRUE~gene2
    ))
    tbl <- table(meta$gene2,meta$Level2)
    row_var <- "gene2"
  }
  
  cell_custom_order <- levels(meta$Level2)
  cell_order <- intersect(cell_custom_order,colnames(tbl))
  obs_mat <- as.data.frame.matrix(tbl)
  obs_mat <- obs_mat[,match(cell_order,colnames(obs_mat))]
  obs_mat <- obs_mat %>% tibble::rownames_to_column(var = "samples")
  
  drop_L1 <- c("Erythrocyte", "Platelet")
  drop_L2_cols <- meta %>% filter(Level1 %in% drop_L1) %>% pull(Level2) %>% unique()
  keep_columns <- setdiff(colnames(obs_mat), drop_L2_cols)
  obs_mat <- obs_mat[, keep_columns]
  
  return(obs_mat)
}

make_scCODA_individual_matrix <- function(
    seurat_obj,
    gene_list,
    fixed_level1_order,
    filter_cell = TRUE
){
  meta_raw <- seurat_obj@meta.data
  meta_raw$Level1 <- factor(meta_raw$Level1, levels = fixed_level1_order)
  level2_order <- meta_raw %>% distinct(Level1, Level2) %>% arrange(Level1, Level2) %>% pull(Level2) %>% unique()
  meta_raw$Level2 <- factor(meta_raw$Level2, levels = level2_order)
  meta_raw <- meta_raw %>% arrange(Level1, Level2)
  
  process_age <- function(df){
    df <- df %>% mutate(age_raw = as.character(age))
    df$age_year <- NA
    idx_m <- str_detect(df$age_raw, regex("m", ignore_case = TRUE))
    df$age_year[idx_m] <- as.numeric(str_remove_all(df$age_raw[idx_m], "[^0-9.]")) / 12
    idx_y <- str_detect(df$age_raw, regex("y", ignore_case = TRUE))
    df$age_year[idx_y] <- as.numeric(str_remove_all(df$age_raw[idx_y], "[^0-9.]"))
    idx_num <- !idx_m & !idx_y & str_detect(df$age_raw, "^[0-9.]+$")
    df$age_year[idx_num] <- as.numeric(df$age_raw[idx_num])
    df$age_group <- cut(
      df$age_year,
      breaks = c(-Inf, 3, 16, 30, 50, 69, Inf),
      labels = c("<3", "3-16", "17-30", "31-50", "51-69", ">69"),
      right = TRUE
    )
    return(df)
  }
  meta_ind <- process_age(meta_raw)
  valid_sample_id <- meta_ind %>% filter(!is.na(sex), !is.na(age_year)) %>% pull(samples) %>% unique()
  meta_use <- meta_ind[meta_ind$samples %in% valid_sample_id, ]
  message(paste0("individual valid samples: ", length(valid_sample_id)))
  tbl <- table(meta_use$samples, meta_use$Level2)
  cell_custom_order <- levels(meta_use$Level2)
  cell_order <- intersect(cell_custom_order, colnames(tbl))
  count_mat <- as.data.frame.matrix(tbl)
  count_mat <- count_mat[, match(cell_order, colnames(tbl))]
  count_mat <- count_mat %>% tibble::rownames_to_column(var = "samples")
  if(filter_cell){
    drop_L1 <- c("Erythrocyte", "Platelet")
    drop_L2_cols <- meta_raw %>% filter(Level1 %in% drop_L1) %>% pull(Level2) %>% unique()
    keep_columns <- setdiff(colnames(count_mat), drop_L2_cols)
    count_mat <- count_mat[, keep_columns]
  }
  cov_info <- meta_use %>% distinct(samples, sex, age_group)
  out_df <- left_join(count_mat, cov_info, by = "samples")
  return(out_df)
}
condition_count<-make_scCODA_count_matrix(pbmc1,"condition",gene_list ,fixed_level1_order)
dim(condition_count);colnames(condition_count)
gene_count<-make_scCODA_count_matrix(pbmc1,"gene",gene_list ,fixed_level1_order)
dim(gene_count);colnames(gene_count)
individual_count<-make_scCODA_individual_matrix(pbmc1,gene_list ,fixed_level1_order)
dim(individual_count);colnames(individual_count)

write.csv(condition_count,file="condition_count.csv",row.names = F)
write.csv(gene_count,file="gene_count.csv",row.names = F)
write.csv(individual_count,file="individual_count.csv",row.names = F)

######2 scCODA results for bubble plot
###(1) continouse
#' Generate scCODA Significant Bubble Plot
#' @title plot_scCODA_bubble
#' @description Visualize scCODA compositional test results with bubble chart.
#' Filter criteria: retain rows with Significant == TRUE AND absolute log2FC >= 2.
#' X axis shows simplified Level2 cell subtype prefix; Y axis shows comparison group label.
#' Point color mapped to raw log2FC: positive value means higher relative abundance in disease group vs HC control.
#' Point size mapped to Inclusion probability (continuous range 0 to 1).
#' Cell subtypes are sorted by predefined fixed Level1 immune lineage order for consistent layout.
#' If no records pass filtering, output blank plot with warning text.
#' @param input_dir Character, directory path storing scCODA output CSV file
#' @param csv_file Character, filename of scCODA result table
#' @returnType ggplot
#' @return ggplot bubble plot object
#' @author MSY
plot_scCODA_bubble <- function(input_dir, csv_file) {
  library(tidyverse)
  
  df_raw <- read_csv(paste0(input_dir, csv_file), show_col_types = FALSE)
  
  # Filter cells with statistical significance and large fold change
  df_filter <- df_raw %>%
    filter(Significant == TRUE, abs(log2FC) >= 2)
  
  # Return empty plot when no qualified records exist
  if (nrow(df_filter) == 0) {
    message("Warning: No cells pass filter (Significant=T & abs(log2FC)>=2)")
    return(ggplot() + annotate("text", x = 1, y = 1, label = "No significant log2FC ≥ 2"))
  }
  
  # Simplify cell subtype label and assign Level1 lineage for fixed ordering
  df_clean <- df_filter %>%
    mutate(
      CellShort = str_extract(`Cell Type`, "^[A-Za-z0-9]+"),
      Level1 = case_when(
        str_detect(CellShort, "^DC") ~ "DC",
        str_detect(CellShort, "^Mono") ~ "Monocyte",
        CellShort %in% c("pDC", "LDG", "Basophil") ~ "Granulocyte/pDC",
        str_detect(CellShort, "^B") ~ "B cell",
        str_detect(CellShort, "^Naive CD4|^NnCD4T") ~ "CD4 T cell",
        str_detect(CellShort, "^Naive CD8|^NnCD8T") ~ "CD8 T cell",
        str_detect(CellShort, "^UTC") ~ "UTC",
        str_detect(CellShort, "^NK") ~ "NK cell",
        TRUE ~ "Other"
      ),
      Level1 = factor(Level1, levels = c(
        "DC",
        "Monocyte",
        "Granulocyte/pDC",
        "B cell",
        "CD4 T cell",
        "CD8 T cell",
        "UTC",
        "NK cell",
        "Other"
      ))
    ) %>%
    arrange(Level1, CellShort) %>%
    mutate(CellShort = factor(CellShort, levels = unique(CellShort)))
  
  # Construct bubble plot
  ggplot(df_clean, aes(x = CellShort, y = Group)) +
    geom_point(
      aes(color = log2FC, size = `Inclusion probability`),
      alpha = 0.85
    ) +
    scale_color_gradient2(
      low = "#4A148C",
      mid = "#FFFFFF",
      high = "#E64A19",
      name = "log2(Fold Change)"
    ) +
    scale_size_continuous(
      range = c(2, 9),
      limits = c(0, 1),
      name = "Inclusion Probability"
    ) +
    theme_bw() +
    theme(
      legend.key.height = unit(0.4, "cm"),
      legend.key.width = unit(0.35, "cm"),
      axis.text = element_text(size = 10, color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      legend.title = element_text(size = 10),
      legend.text = element_text(color = "black"),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    labs(x = "Cell Subtype", y = "")
}
plot_scCODA_bubble(
  input_dir = "./",
  csv_file = "scCODA_condition_out.csv")

###(2) discreat
#' Generate scCODA Bubble Plot V2
#' @title plot_scCODA_bubble2
#' @description Visualize filtered scCODA compositional test bubble chart.
#' Filter rule: Significant == TRUE & abs(log2FC) >= adjustable fc_cutoff.
#' X axis: simplified Level2 cell prefix sorted by fixed immune Level1 lineage.
#' Y axis: comparison group order reversed to flip vertical display.
#' Point size mapped to Inclusion probability 3 discrete tiers; color mapped to 6 discrete log2FC magnitude tiers.
#' Vertical dashed lines only at log2FC = -1 / 1; inclusion cutoff horizontal dashed line color #898989.
#' Light gray panel grid, remove x/y axis labels, all text/ticks pure black, no point stroke border.
#' @param input_dir Character, folder path of scCODA csv file
#' @param csv_file Character, scCODA result csv filename
#' @param fc_cutoff Numeric, absolute log2FC filter threshold, default = 0.5
#' @param tier_color_vec Character vector length=6, custom six-tier log2FC color palette, NULL = built-in cold-warm palette
#' @returnType ggplot
#' @author MSY
plot_scCODA_bubble2 <- function(input_dir,
                                csv_file,
                                fc_cutoff = 0.5,
                                tier_color_vec = NULL) {
  library(tidyverse)
  
  # Original screenshot top-to-bottom group sequence
  group_raw_order <- c(
    "KD",
    "JIA",
    "SLE",
    "Uncategorized",
    "Protein_homeostasis",
    "TBK1_IRF3",
    "Cytoskeleton_and_small_GTPase",
    "Negative_regulation_of_IFN-I",
    "Inborn_errors_of_cell_death",
    "Immune_metabolic",
    "Inflammasome_IL-1β",
    "Arachidonic acid metabolism",
    "Endolysosomal_nucleic_acid_sensing",
    "Ca2+_flux-PLC",
    "Osteoclast function",
    "NFKB pathway"
  )
  # Reverse levels to flip Y axis display
  group_levels <- rev(group_raw_order)
  
  # Built-in six discrete tier color: 3 cold negative, 3 warm positive
  default_tier_palette <- c(
    "< -2" = "#283593",
    "-2 ~ -1" = "#5C6BC0",
    "-1 ~ 0" = "#BBDEFB",
    "0 ~ 1" = "#FFCCBC",
    "1 ~ 2" = "#E53935",
    "> 2" = "#B71C1C"
  )
  # 
  # default_tier_palette <- c(
  #   "< -2" = "#003D81",
  #   "-2 ~ -1" = "#0B71AB",
  #   "-1 ~ 0" = "#9CD2ED",
  #   "0 ~ 1" = "#F7BE6F",
  #   "1 ~ 2" = "#E55C3A",
  #   "> 2" = "#8C2522"
  # )
  # default_tier_palette <- c(
  #   "< -2" = "#003D81",
  #   "-2 ~ -1" = "#5C6BC0",
  #   "-1 ~ 0" = "#BBDEFB",
  #   "0 ~ 1" = "#FFCCBC",
  #   "1 ~ 2" = "#E53935",
  #   "> 2" ="#8C2522"
  # )
  
  
  
  #003D81 "#0B71AB" "#9CD2ED" "#F7BE6F" "#E55C3A" #8C2522"
  df_raw <- read_csv(paste0(input_dir, csv_file), show_col_types = FALSE)
  df_filter <- df_raw %>% filter(Significant == TRUE, abs(log2FC) >= fc_cutoff)
  
  # Empty plot template
  if (nrow(df_filter) == 0) {
    message(paste0("Warning: No cells pass filter (Significant=T & abs(log2FC) >= ", fc_cutoff, ")"))
    return(ggplot() +
             annotate("text", x = 1, y = 0.5, label = paste0("No cells meet abs(log2FC) ≥ ", fc_cutoff)) +
             labs(x = "", y = "") +
             theme_bw() +
             theme(
               panel.grid = element_line(color = "gray90", linewidth = 0.2),
               axis.text = element_text(color = "black"),
               axis.ticks = element_line(color = "black"),
               legend.position = "none"
             ))
  }
  
  df_clean <- df_filter %>%
    mutate(
      CellShort = str_extract(`Cell Type`, "^[A-Za-z0-9]+"),
      Level1 = case_when(
        str_detect(CellShort, "^DC") ~ "DC",
        str_detect(CellShort, "^Mono") ~ "Monocyte",
        CellShort %in% c("pDC", "LDG", "Basophil") ~ "Granulocyte/pDC",
        str_detect(CellShort, "^B") ~ "B cell",
        str_detect(CellShort, "^Naive CD4|^NnCD4T") ~ "CD4 T cell",
        str_detect(CellShort, "^Naive CD8|^NnCD8T") ~ "CD8 T cell",
        str_detect(CellShort, "^UTC") ~ "UTC",
        str_detect(CellShort, "^NK") ~ "NK cell",
        TRUE ~ "Other"
      ),
      Level1 = factor(Level1, levels = c("DC","Monocyte","Granulocyte/pDC","B cell","CD4 T cell","CD8 T cell","UTC","NK cell","Other")),
      # Inclusion probability three size tiers
      IncTier = case_when(
        `Inclusion probability` < 0.5 ~ "Low (<0.5)",
        `Inclusion probability` >= 0.5 & `Inclusion probability` < 0.8 ~ "Medium (0.5~0.8)",
        `Inclusion probability` >= 0.8 ~ "High (≥0.8)"
      ),
      # Six log2FC magnitude tiers for discrete color
      FcTier = case_when(
        log2FC < -2 ~ "< -2",
        log2FC >= -2 & log2FC < -1 ~ "-2 ~ -1",
        log2FC >= -1 & log2FC < 0 ~ "-1 ~ 0",
        log2FC >= 0 & log2FC < 1 ~ "0 ~ 1",
        log2FC >= 1 & log2FC < 2 ~ "1 ~ 2",
        log2FC >= 2 ~ "> 2"
      ),
      Group = factor(Group, levels = group_levels)
    ) %>%
    arrange(Level1, CellShort) %>%
    mutate(CellShort = factor(CellShort, levels = unique(CellShort)))
  
  ggplot(df_clean, aes(x = CellShort, y = Group)) +
    # Vertical dashed reference lines only at -1 / 1
    geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "gray60", linewidth = 0.3) +
    # Inclusion threshold horizontal dashed line #898989
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "#898989", linewidth = 0.4) +
    geom_point(aes(color = FcTier, size = IncTier), alpha = 0.85, stroke = 0) +
    scale_color_manual(values = if(!is.null(tier_color_vec)) setNames(tier_color_vec, names(default_tier_palette)) else default_tier_palette, name = "log2FC Tier") +
    scale_size_manual(values = c("Low (<0.5)" = 2, "Medium (0.5~0.8)" = 4, "High (≥0.8)" = 6), name = "Inclusion Probability") +
    theme_bw() +
    theme(
      panel.grid = element_line(color = "gray90", linewidth = 0.2),
      legend.key.height = unit(0.4, "cm"),
      legend.key.width = unit(0.35, "cm"),
      axis.text = element_text(size = 10, color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_text(size = 10),
      legend.text = element_text(color = "black"),
      plot.margin = margin(5,5,5,5)
    ) +
    labs(x = "", y = "")
}


p_05 <- plot_scCODA_bubble2(
  input_dir = "./",
  csv_file = "scCODA_condition_out.csv",
  fc_cutoff = 0.5
)
p_05
p_10 <- plot_scCODA_bubble2(
  input_dir = "./",
  csv_file = "scCODA_condition_out.csv",
  fc_cutoff = 1
)


#######3 scCODA for scatter 
#' Generate Scatter Plot List for Each scCODA Comparison Group
#' @title plot_scCODA_scatter_allgroups
#' @description Create individual scatter plot for every unique Group in scCODA result table.
#' Only retain cells with Significant == TRUE for each subplot.
#' X axis: log2FC; Y axis: Inclusion probability.
#' Vertical dashed lines at x = -1, 1 for fold change reference.
#' Fixed horizontal threshold dashed line (#898989) for Inclusion probability, adjustable via inc_cutoff parameter.
#' Points colored by predefined Level1 immune lineage palette, all points same size.
#' Points with abs(log2FC) > 1 AND Inclusion probability >= inc_cutoff get black stroke & text labels.
#' All plots stored in named list for manual selection & export.
#' Hide color legend, plot title only group name and centered. All axis text/ticks set to black.
#' Remove x & y axis labels to save layout space when stitching.
#' @param input_dir Character, path directory of scCODA csv result
#' @param csv_file Character, scCODA output table filename
#' @param inc_cutoff Numeric, fixed threshold of inclusion probability (0~1), default = 0.6
#' @returnType List
#' @return Named list of ggplot objects, names equal to each comparison Group
#' @author MSY
plot_scCODA_scatter_allgroups <- function(input_dir,
                                          csv_file,
                                          inc_cutoff = 0.8) {
  library(tidyverse)
  
  # Fixed immune lineage color palette provided
  my_colors <- c(
    "Naive CD4 T"="#847AB3",
    "Non-naive CD4 T"="#C8BBD5",
    "Naive CD8 T"="#F6E1EE",
    "Non-naive CD8 T"="#0B71AB",
    "NK"="#96AF95",
    "B cell"="#003D81",
    "Monocyte"="#F49D5C",
    "LDG"="#C95968",
    "Plasma"="#B3D1E7",
    "Erythrocyte"="#FFEFC1",
    "Platelet"="#B3928B",
    "UTC"="#86C7B4",
    "DC"="#F0E442",
    "pDC"="#B07A99",
    "Basophil"="#8C2522"
  )
  
  df_raw <- read_csv(paste0(input_dir, csv_file), show_col_types = FALSE)
  
  # Preprocess full dataset: simplify cell name & assign Level1 lineage & color group
  df_all <- df_raw %>%
    mutate(
      CellShort = str_extract(`Cell Type`, "^[A-Za-z0-9]+"),
      Level1 = case_when(
        str_detect(CellShort, "^DC") ~ "DC",
        str_detect(CellShort, "^Mono") ~ "Monocyte",
        CellShort == "LDG" ~ "LDG",
        CellShort == "Basophil" ~ "Basophil",
        CellShort == "pDC" ~ "pDC",
        str_detect(CellShort, "^B") ~ "B cell",
        str_detect(CellShort, "^Naive CD4") ~ "Naive CD4 T",
        str_detect(CellShort, "^NnCD4T") ~ "Non-naive CD4 T",
        str_detect(CellShort, "^Naive CD8") ~ "Naive CD8 T",
        str_detect(CellShort, "^NnCD8T") ~ "Non-naive CD8 T",
        str_detect(CellShort, "^UTC") ~ "UTC",
        str_detect(CellShort, "^NK") ~ "NK",
        TRUE ~ "Other"
      )
    )
  
  unique_groups <- unique(df_all$Group)
  plot_list <- list()
  
  # Loop over every single comparison group
  for (grp in unique_groups) {
    df_sub <- df_all %>% filter(Group == grp, Significant == TRUE)
    
    if (nrow(df_sub) == 0) {
      p_empty <- ggplot() +
        annotate("text", x = 0, y = 0.5, label = paste0(grp, ": No significant cells")) +
        labs(title = grp, x = "", y = "") +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = 0.5, color = "black"),
          axis.text = element_text(color = "black"),
          axis.ticks = element_line(color = "black"),
          legend.position = "none"
        )
      plot_list[[grp]] <- p_empty
      next
    }
    
    # Flag points that need black stroke & text label
    df_sub <- df_sub %>%
      mutate(
        highlight = case_when(
          abs(log2FC) > 1 & `Inclusion probability` >= inc_cutoff ~ TRUE,
          TRUE ~ FALSE
        )
      )
    
    p <- ggplot(df_sub, aes(x = log2FC, y = `Inclusion probability`)) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray60", linewidth = 0.3) +
      geom_hline(yintercept = inc_cutoff, linetype = "dashed", color = "#898989", linewidth = 0.4) +
      # All points fixed same size=3, only stroke differs
      geom_point(aes(color = Level1, stroke = ifelse(highlight, 1.2, 0)), alpha = 0.8, size = 3) +
      scale_color_manual(values = my_colors, name = "Lineage") +
      geom_text(
        data = filter(df_sub, highlight == TRUE),
        aes(label = CellShort),
        size = 2.7,
        hjust = -0.15,
        color = "black"
      ) +
      # Remove x/y axis labels
      labs(x = "", y = "", title = grp) +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        plot.title = element_text(size = 11, hjust = 0.5, color = "black"),
        axis.text = element_text(size = 9, color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none",
        plot.margin = margin(5, 5, 5, 5)
      )
    
    plot_list[[grp]] <- p
  }
  
  return(plot_list)
}

#' Combine multiple scCODA scatter plots into one grid figure
#' @title combine_scCODA_scatter_grid
#' @description Arrange all subplots from plot list into grid layout.
#' All subplots already removed x/y axis titles internally.
#' @param plot_list Named list of ggplot objects output from plot_scCODA_scatter_allgroups
#' @param ncol Integer, number of columns in grid layout
#' @returnType gtable
#' @return Combined grid plot object
#' @author MSY
combine_scCODA_scatter_grid <- function(plot_list, ncol = 4) {
  library(patchwork)
  wrap_plots(plot_list, ncol = ncol)
}

plot_collection <- plot_scCODA_scatter_allgroups(
  input_dir = "./",
  csv_file = "scCODA_condition_out.csv",
  inc_cutoff = 0.6)
all_grid <- combine_scCODA_scatter_grid(plot_collection, ncol = 4)
all_grid




