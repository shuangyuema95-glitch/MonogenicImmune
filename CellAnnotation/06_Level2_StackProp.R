library(Seurat)
library(tidyverse)
library(scales)

plot_subtype_stacked <- function(
    seu_obj,
    type = c("sample", "gene", "condition"),
    target_cell,
    gene_list,
    meta_col_level1 = "Level1",
    meta_col_level2 = "Level2",
    meta_col_sample = "samples",
    meta_col_gene = "gene",
    condition_order = c("HC", "KD", "JIA","SLE","Uncategorized","Protein_homeostasis",
                        "TBK1_IRF3","Cytoskeleton_and_small_GTPase","Negative_regulation_of_IFN-I","Inborn_errors_of_cell_death",
                        "Immune_metabolic","Inflammasome_IL-1β","Arachidonic acid metabolism",
                        "Endolysosomal_nucleic_acid_sensing",
                        "Ca2+_flux-PLC","Osteoclast function","NFKB pathway"),
    show_pct_label = FALSE
){
  type <- match.arg(type)
  meta <- seu_obj@meta.data %>% as_tibble()
  
  meta <- meta %>% filter(.data[[meta_col_level1]] == target_cell)
  if(nrow(meta) == 0) stop("No cells for target_cell in Level1")
  
  gene_map <- stack(gene_list)
  colnames(gene_map) <- c("gene", "pathway")
  
  if(type == "sample"){
    meta <- meta %>% mutate(group = .data[[meta_col_sample]])
    row_var <- "group"
    x_levels <- sort(unique(meta$group))
  } else if(type == "condition"){
    meta$condition0 <- gene_map$pathway[match(meta[[meta_col_gene]], gene_map$gene)]
    meta <- meta %>% mutate(
      group = case_when(
        str_detect(.data[[meta_col_sample]], "SLE") ~ "SLE",
        str_detect(.data[[meta_col_sample]], "JIA") ~ "JIA",
        str_detect(.data[[meta_col_sample]], "KD") ~ "KD",
        str_detect(.data[[meta_col_sample]], "HC") ~ "HC",
        TRUE ~ condition0
      )
    )
    row_var <- "group"
    x_levels <- intersect(condition_order, unique(meta$group))
  } else if(type == "gene"){
    meta$pathway <- gene_map$pathway[match(meta[[meta_col_gene]], gene_map$gene)]
    meta <- meta %>% mutate(
      group = case_when(
        str_detect(.data[[meta_col_sample]], "SLE") ~ "SLE",
        str_detect(.data[[meta_col_sample]], "JIA") ~ "JIA",
        str_detect(.data[[meta_col_sample]], "KD") ~ "KD",
        str_detect(.data[[meta_col_sample]], "HC") ~ "HC",
        TRUE ~ .data[[meta_col_gene]]
      )
    )
    row_var <- "group"
    pathway_vec <- condition_order[! condition_order %in% c("SLE","JIA","KD","HC")]
    gene_order <- c()
    for(pw in pathway_vec){
      if(pw %in% names(gene_list)){
        gene_order <- c(gene_order, gene_list[[pw]])
      }
    }
    # gene模式：先HC,KD,JIA,SLE，再接各个基因
    x_levels <- c(c("HC","KD","JIA","SLE"), gene_order)
    x_levels <- intersect(x_levels, unique(meta$group))
  }
  
  count_df <- meta %>%
    filter(!is.na(group)) %>%
    count(.data[[row_var]], .data[[meta_col_level2]], name = "n_sub") %>%
    group_by(.data[[row_var]]) %>%
    mutate(n_total = sum(n_sub), frac = n_sub / n_total) %>%
    ungroup()
  
  lv2_order_df <- count_df %>% distinct(.data[[meta_col_level2]]) %>%
    mutate(num_tag = str_extract(.data[[meta_col_level2]], "(?<=0|_)\\d+") %>% as.integer()) %>%
    arrange(num_tag)
  lv2_levels <- lv2_order_df[[meta_col_level2]]
  
  count_df[[meta_col_level2]] <- factor(count_df[[meta_col_level2]], levels = lv2_levels)
  count_df[[row_var]] <- factor(count_df[[row_var]], levels = x_levels)
  
  n_subtypes <- length(lv2_levels)
  set.seed(123)
  subtype_cols <- hue_pal(h = c(0,360), l = 65, c = 85)(n_subtypes)
  names(subtype_cols) <- lv2_levels
  
  p <- ggplot(count_df, aes(x = .data[[row_var]], y = frac, fill = .data[[meta_col_level2]])) +
    geom_col(position = position_stack(reverse = FALSE), width = 0.75) +
    scale_fill_manual(values = subtype_cols) +
    scale_y_continuous(expand = c(0, 0), labels = percent_format()) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
      axis.text.y = element_text(colour = "black"),
      legend.title = element_blank(),
      panel.grid = element_blank()
    ) +
    labs(x = "", y = str_glue("Fraction within {target_cell}"))
  
  if(isTRUE(show_pct_label)){
    p <- p + geom_text(
      aes(label = ifelse(frac > 0.01, percent(frac, accuracy = 1), "")),
      position = position_stack(vjust = 0.5),
      colour = "white", size = 3.2
    )
  }
  
  return(list(plot = p, prop_table = count_df))
}

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
  "Protein_homeostasis" = c("UBA1", "PSMD12")
)

# example, show percentage label
res_gene <- plot_subtype_stacked(
  seu_obj = pbmc1,
  type = "gene",
  target_cell = "Non-naive CD8 T",
  gene_list = gene_list,
  show_pct_label = TRUE
)
res_gene$plot

res_condition <- plot_subtype_stacked(
  seu_obj = pbmc1,
  type = "condition",
  target_cell = "Non-naive CD8 T",
  gene_list = gene_list,
  show_pct_label = FALSE
)
res_condition$plot

res_sample <- plot_subtype_stacked(
  seu_obj = pbmc1,
  type = "sample",
  target_cell = "Non-naive CD8 T",
  gene_list = gene_list,
  show_pct_label = FALSE
)
res_sample$plot