# ============================================================
# LMEM: t_stat ~ disease + chemistry + (1|studyID), 218 factors x 16 diseases vs HC.
# Heatmaps: continuous + 8-level discrete, each with cells/immunecate/pathway modes.
# Output: lmem_results_norare.csv + lmem_results_annotated.xlsx/txt.
# beta = disease-vs-HC effect size in ULM t-stat units (positive = up in disease).
# ============================================================


library(dplyr)
library(readr)
library(lme4)
library(lmerTest)
library(stringr)
library(tidyr)
library(data.table)

setwd("E:\\AID cohort\\code\\NMF\\Spectra")

# ============================================================
# LMEM without rare-disease merging (17 disease categories)
# Formula: t_stat ~ disease + chemistry + (1 | studyID)
# ============================================================
library(dplyr); library(readr); library(data.table)
library(lme4); library(lmerTest); library(stringr); library(tidyr)
library(pheatmap)
setwd("E:\\AID cohort\\code\\NMF\\Spectra")

######STEP1: load ULM result + metadata
ulm_res <- read.csv("lam0001_dict242_second_ulm_result.csv", stringsAsFactors = FALSE)
meta <- as.data.frame(fread("E:\\AID cohort\\code\\metainfo.txt"))

gene_list <- list(
  "NFKB pathway" = c("TNFAIP3","IRAK2","NOD2","RELA"),
  "Uncategorized" = c("ADA2","ELF4","CSF3R","STAT4"),
  "Osteoclast function" = c("ACP5","OGFRL1"),
  "Ca2+_flux-PLC" = c("PLCG2","PLCG1"),
  "Endolysosomal_nucleic_acid_sensing" = c("UNC93B1","TLR8","TLR7","TLR1","PLD4","TRAF3"),
  "Arachidonic acid metabolism" = c("TBXAS1"),
  "Inflammasome_IL-1β" = c("IL1R1","NLRC4","NLRP3","PSTPIP1","CDC42","LPIN2"),
  "Immune_metabolic" = c("LACC1","SLC7A7"),
  "Inborn_errors_of_cell_death" = c("RIPK1","TNFRSF1A","OTULIN","RNF31"),
  "Negative_regulation_of_IFN-I" = c("ISG15","USP18","SOCS1"),
  "Cytoskeleton_and_small_GTPase" = c("GNAI2","KRAS"),
  "TBK1_IRF3" = c("STING","TREX1","IFIH1","COPA"),
  "Protein_homeostasis" = c("UBA1","PSMD12"))
gene_map <- stack(gene_list); colnames(gene_map) <- c("gene","condition0")

meta$samples <- meta$dataset
meta$condition0 <- gene_map$condition0[match(meta$gene, gene_map$gene)]
meta <- meta %>% mutate(condition = case_when(
  str_detect(samples,"SLE") ~ "SLE",
  str_detect(samples,"JIA") ~ "JIA",
  str_detect(samples,"KD") ~ "KD",
  str_detect(samples,"HC") ~ "HC",
  TRUE ~ condition0))

# NO rare merging: keep all 13 monogenic mechanisms + 3 polygenic + HC
meta$disease <- factor(meta$condition, levels = c(
  "HC","SLE","Inborn_errors_of_cell_death","TBK1_IRF3",
  "Endolysosomal_nucleic_acid_sensing","Inflammasome_IL-1β","NFKB pathway",
  "JIA","KD","Uncategorized","Negative_regulation_of_IFN-I",
  "Ca2+_flux-PLC","Osteoclast function","Arachidonic acid metabolism",
  "Cytoskeleton_and_small_GTPase","Immune_metabolic","Protein_homeostasis"))
table(meta$disease)


######STEP2: merge + factorize
ulm_res <- ulm_res %>% left_join(meta, by = "samples")
ulm_res$chemistry <- factor(ulm_res$chemistry)
ulm_res$studyID <- factor(ulm_res$resource)

######STEP3: celltype-specific filtering
extract_ct <- function(x) str_match(x, "^Factor_\\d+[ _]([^_ ]+)")[, 2]
ulm_res$factor_ct <- sapply(ulm_res$factor, extract_ct)

ulm_res$level2_ct <- case_when(
  ulm_res$Level2 == "Basophil" ~ "Baso",
  str_detect(ulm_res$Level2, "^B") ~ "B",
  str_detect(ulm_res$Level2, "^DC") ~ "DC",
  ulm_res$Level2 == "LDG" ~ "LDG",
  str_detect(ulm_res$Level2, "^Mono") ~ "Mono",
  str_detect(ulm_res$Level2, "^NK") ~ "NK",
  ulm_res$Level2 == "Naive CD4 T" ~ "NnCD4T",
  ulm_res$Level2 == "Naive CD8 T" ~ "NnCD8T",
  str_detect(ulm_res$Level2, "^NnCD4T") ~ "NonNCD4T",
  str_detect(ulm_res$Level2, "^NnCD8T") ~ "NonNCD8T",
  str_detect(ulm_res$Level2, "^Plasma") ~ "Plasma",
  ulm_res$Level2 == "Platelet" ~ "Platelet",
  str_detect(ulm_res$Level2, "^UTC") ~ "UTC",
  ulm_res$Level2 == "pDC" ~ "pDC",
  TRUE ~ "OTHER")

ulm_filtered <- ulm_res %>% filter(factor_ct == "all" | factor_ct == level2_ct)
cat("Filtered:", dim(ulm_filtered), "factors:", length(unique(ulm_filtered$factor)), "\n")

######: LMEM per factor
disease_levels <- setdiff(levels(ulm_filtered$disease), "HC")
factors <- unique(ulm_filtered$factor)
lmem_res <- data.frame()

for (i in seq_along(factors)) {
  f <- factors[i]
  df_sub <- ulm_filtered %>% filter(factor == !!f)
  ct <- unique(df_sub$factor_ct)
  n_hc <- sum(df_sub$disease == "HC", na.rm = TRUE)
  
  if (n_hc < 3) {
    for (d in disease_levels)
      lmem_res <- rbind(lmem_res, data.frame(
        factor=f, celltype=ct, disease=d, n_hc=n_hc, n_dis=NA,
        beta=NA, se=NA, t_val=NA, p_value=NA, converged=FALSE, note="HC<3"))
    next
  }
  
  fit <- tryCatch({
    lmer(t_stat ~ disease + chemistry + (1 | studyID), data = df_sub,
         control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
  }, error = function(e) NULL)
  
  if (is.null(fit)) {
    for (d in disease_levels)
      lmem_res <- rbind(lmem_res, data.frame(
        factor=f, celltype=ct, disease=d, n_hc=n_hc,
        n_dis=sum(df_sub$disease==d, na.rm=TRUE),
        beta=NA, se=NA, t_val=NA, p_value=NA, converged=FALSE, note="model error"))
    next
  }
  
  coef_tab <- as.data.frame(summary(fit)$coefficients)
  for (r in grep("^disease", rownames(coef_tab))) {
    d <- sub("^disease", "", rownames(coef_tab)[r])
    n_dis <- sum(df_sub$disease == d, na.rm = TRUE)
    lmem_res <- rbind(lmem_res, data.frame(
      factor=f, celltype=ct, disease=d, n_hc=n_hc, n_dis=n_dis,
      beta=coef_tab[r,"Estimate"], se=coef_tab[r,"Std. Error"],
      t_val=coef_tab[r,"t value"], p_value=coef_tab[r,"Pr(>|t|)"],
      converged=TRUE, note=ifelse(n_dis<3, "dis<3", "")))
  }
  if (i %% 50 == 0) cat(sprintf("  [%d/%d] %s\n", i, length(factors), f))
}
cat(sprintf("  [%d/%d] done\n", length(factors), length(factors)))

######: FDR + save
lmem_res$fdr_global <- p.adjust(lmem_res$p_value, method = "fdr")
lmem_res <- lmem_res %>% group_by(celltype) %>%
  mutate(fdr_celltype = p.adjust(p_value, method = "fdr")) %>% ungroup()
lmem_res <- lmem_res %>% arrange(fdr_global)
write.csv(lmem_res, "lam0001_dict242_lmem_results_norare.csv", row.names = FALSE)

cat("Total pairs:", nrow(lmem_res), "\n")
cat("Converged:", sum(lmem_res$converged), "\n")
cat("Sig (fdr_global<0.05):", sum(lmem_res$fdr_global < 0.05, na.rm=TRUE), "\n")
cat("NA beta:", sum(is.na(lmem_res$beta)), "\n")
lmem_res %>% filter(fdr_global < 0.05) %>% group_by(disease) %>%
  summarise(n_sig=n(), n_up=sum(beta>0), n_down=sum(beta<0)) %>% print()

###### STEP6: heatmap function (two modes)
###(1)continuous color value for beta values produced by LMEM model
plot_lmem_heatmap <- function(lmem_res,
                              mode = c("cells", "immunecate", "pathway"),
                              factor_list = NULL,
                              strip_celltype_prefix = FALSE,
                              cellwidth = NA,
                              cellheight = 12,
                              heatmap_cols = c("#5E4B8B", "#C8BBD5", "white", "#E8C4A0", "#8C2522"),
                              immune_rds = "E:\\AID cohort\\code\\NMF\\Spectra\\0830_PieImmuneCate.RDS",
                              na_text_size = 4,
                              title = NULL) {
  mode <- match.arg(mode)
  if (is.null(title)) title <- sprintf("%s_groupby%s (LMEM beta, disease vs HC)",
                                       ifelse(is.null(factor_list), "whole", "partial"), mode)
  library(pheatmap); library(dplyr); library(tidyr); library(stringr); library(gtable); library(grid)
  
  # matrices
  beta_mat <- lmem_res %>% select(disease, factor, beta) %>%
    pivot_wider(names_from = factor, values_from = beta) %>% column_to_rownames("disease") %>% as.matrix()
  sig_mat <- lmem_res %>% select(disease, factor, fdr_global) %>%
    mutate(sig = ifelse(fdr_global < 0.05, "*", "")) %>% select(-fdr_global) %>%
    pivot_wider(names_from = factor, values_from = sig) %>% column_to_rownames("disease") %>% as.matrix()
  
  # palettes
  ct_order <- c("DC","Mono","pDC","LDG","Baso","B","Plasma",
                "NnCD4T","NonNCD4T","NnCD8T","NonNCD8T","UTC","NK","Platelet","all")
  celltype_colors <- c("DC"="#F0E442","Mono"="#F49D5C","pDC"="#B07A99","LDG"="#C95968",
                       "Baso"="#8C2522","B"="#003D81","Plasma"="#B3D1E7","NnCD4T"="#847AB3","NonNCD4T"="#C8BBD5",
                       "NnCD8T"="#F6E1EE","NonNCD8T"="#0B71AB","UTC"="#86C7B4","NK"="#96AF95","Platelet"="#B3928B","all"="#502680")
  Immune_cateOrder <- c("Cell function and homeostasis","Inflammatory regulation","Adaptive immune regulation",
                        "Innate immune sensing","Amino acid metabolism","Energy metabolism","Lipid metabolism",
                        "Small molecule & nucleotide metabolism","Intracellular signal transduction")
  imm_colors <- setNames(c("#2F91A7","#7DB6C1","#C1DCE2","#8E9FC9","#D5B3CD","#AA94C1","#E9B48F","#E1CF84","#E4935D"), Immune_cateOrder)
  sp_colors <- c("cell-type-specific"="#7777BD", "Non cell-type-specific"="#CF7D83")
  group_colors <- c("polygenic patient"="#B18749", "monogenic patient"="#8C2522")
  
  disease_colors <- c("KD"="#1976D2","JIA"="#7F7FD5","SLE"="#998EC3","Uncategorized"="#B07AA1",
                      "Protein_homeostasis"="#DC9497","TBK1_IRF3"="#9D7660","Cytoskeleton_and_small_GTPase"="#BAB0AC",
                      "Negative_regulation_of_IFN-I"="#76B7B2","Inborn_errors_of_cell_death"="#F28E2B",
                      "Immune_metabolic"="#EDC948",
                      "Inflammasome_IL-1β"= "#8A9FD1","Arachidonic acid metabolism"="#61BADA",
                      "Endolysosomal_nucleic_acid_sensing"="#D37295",
                      "Ca2+_flux-PLC"="#A0CBE8","Osteoclast function"="#FABFD2","NFKB pathway"="#FFBE7D")
  
  # col_anno
  factor_ct_vec <- setNames(lmem_res$celltype, lmem_res$factor)
  all_factors <- colnames(beta_mat)
  col_anno <- data.frame(factor_id = all_factors, celltype = factor_ct_vec[all_factors],
                         specificity = ifelse(factor_ct_vec[all_factors] == "all", "Non cell-type-specific", "cell-type-specific"),
                         name = sub("^Factor_\\d+[ _]", "", all_factors), row.names = all_factors, stringsAsFactors = FALSE)
  col_anno$pure_pathway <- sub("^[^_]+_", "", col_anno$name)
  ImmuneCate <- readRDS(immune_rds)
  ct_prefix_regex <- paste0("^(", paste(rev(ct_order), collapse="|"), ")_")
  col_anno$immune_category <- ImmuneCate$new_category_en[match(sub(ct_prefix_regex, "", col_anno$name), ImmuneCate$name)]
  
  # partial filter
  if (!is.null(factor_list)) {
    sel <- factor_list[factor_list != ""]
    matched <- col_anno$factor_id %in% sel | col_anno$name %in% sel
    cat(sprintf("[partial] %d / %d factors matched\n", sum(matched), length(sel)))
    col_anno <- col_anno[matched, ]; beta_mat <- beta_mat[, rownames(col_anno)]; sig_mat <- sig_mat[, rownames(col_anno)]
  }
  
  # sort + gaps by mode
  if (mode == "cells") {
    col_anno <- col_anno[order(factor(col_anno$celltype, levels = ct_order)), ]
    ct_counts <- table(factor(col_anno$celltype, levels = ct_order)); ct_counts <- ct_counts[ct_counts > 0]
    gaps_col <- cumsum(ct_counts)[-length(ct_counts)]
  } else if (mode == "immunecate") {
    col_anno <- col_anno[order(factor(col_anno$immune_category, levels = Immune_cateOrder), factor(col_anno$celltype, levels = ct_order)), ]
    imm_counts <- table(factor(col_anno$immune_category, levels = Immune_cateOrder)); imm_counts <- imm_counts[imm_counts > 0]
    gaps_col <- cumsum(imm_counts)[-length(imm_counts)]
  } else {
    pw_order <- names(sort(table(col_anno$pure_pathway), decreasing = TRUE))
    col_anno <- col_anno[order(factor(col_anno$pure_pathway, levels = pw_order), factor(col_anno$celltype, levels = ct_order)), ]
    pw_counts <- table(factor(col_anno$pure_pathway, levels = pw_order)); pw_counts <- pw_counts[pw_counts > 0]
    gaps_col <- cumsum(pw_counts)[-length(pw_counts)]
  }
  beta_mat <- beta_mat[, rownames(col_anno)]; sig_mat <- sig_mat[, rownames(col_anno)]
  colnames(beta_mat) <- col_anno$name; colnames(sig_mat) <- col_anno$name; rownames(col_anno) <- col_anno$name
  col_anno_plot <- col_anno[, c("celltype", "specificity", "immune_category")]
  
  # column labels
  if (mode == "pathway") {
    labels_col <- rep("", nrow(col_anno))
  } else if (strip_celltype_prefix) {
    labels_col <- sub("^[^_]+_", "", colnames(beta_mat))
  } else { labels_col <- colnames(beta_mat) }
  
  # row order + annotation (group + disease as color bars, no text)
  condition_order <- c("KD","JIA","SLE","Uncategorized","Protein_homeostasis","TBK1_IRF3",
                       "Cytoskeleton_and_small_GTPase","Negative_regulation_of_IFN-I","Inborn_errors_of_cell_death",
                       "Immune_metabolic","Inflammasome_IL-1β","Arachidonic acid metabolism","Endolysosomal_nucleic_acid_sensing",
                       "Ca2+_flux-PLC","Osteoclast function","NFKB pathway")
  beta_mat <- beta_mat[condition_order, ]; sig_mat <- sig_mat[condition_order, ]
  disease_group <- data.frame(group = c(rep("polygenic patient", 3), rep("monogenic patient", 13)),
                              disease = condition_order, row.names = condition_order)
  
  # annotation colors (keep only existing categories)
  actual_ct <- unique(col_anno$celltype); actual_imm <- unique(col_anno$immune_category); actual_sp <- unique(col_anno$specificity)
  ann_colors <- list(group = group_colors, disease = disease_colors[condition_order],
                     celltype = celltype_colors[names(celltype_colors) %in% actual_ct],
                     specificity = sp_colors[names(sp_colors) %in% actual_sp],
                     immune_category = imm_colors[names(imm_colors) %in% actual_imm])
  
  # heatmap
  bk <- seq(-max(abs(beta_mat), na.rm=TRUE), max(abs(beta_mat), na.rm=TRUE), length.out=100)
  cols <- colorRampPalette(heatmap_cols)(100)
  n_col <- ncol(beta_mat); fontsize_col <- ifelse(n_col > 150, 5, ifelse(n_col > 100, 6, 7))
  p <- pheatmap(beta_mat, color = cols, breaks = bk, cluster_rows = FALSE, cluster_cols = FALSE,
                gaps_row = 3, gaps_col = gaps_col, annotation_col = col_anno_plot, annotation_row = disease_group,
                annotation_colors = ann_colors, annotation_names_row = FALSE, show_rownames = FALSE,
                labels_col = labels_col, display_numbers = sig_mat, number_color = "black", fontsize_number = 7,
                border_color = "white", main = title, angle_col = 90, fontsize_row = 10, fontsize_col = fontsize_col,
                legend_breaks = c(-max(abs(beta_mat),na.rm=TRUE), 0, max(abs(beta_mat),na.rm=TRUE)),
                legend_labels = c("Down", "0", "Up"), cellwidth = cellwidth, cellheight = cellheight, silent = TRUE)
  
  # shrink annotation legend color blocks
  ann_leg_idx <- grep("annotation_legend", p$gtable$layout$name)
  if (length(ann_leg_idx) > 0) {
    for (cl in unique(p$gtable$layout$l[ann_leg_idx])) p$gtable$widths[cl] <- unit(0.85, "cm")
    for (i in ann_leg_idx) {
      g <- editGrob(p$gtable$grobs[[i]], gPath("rect"), grep = TRUE, global = TRUE,
                    height = unit(0.32, "cm"), width = unit(0.32, "cm"))
      p$gtable$grobs[[i]] <- g
    }
  }
  
  # shrink NA text (stars stay at 7)
  shrink_na <- function(g) {
    if (is.grob(g)) {
      if (inherits(g, "text") && !is.null(g$label) && ("NA" %in% g$label)) g$gp$fontsize <- na_text_size
      if (!is.null(g$children)) g$children <- lapply(g$children, shrink_na)
    }
    return(g)
  }
  p$gtable$grobs <- lapply(p$gtable$grobs, shrink_na)
  
  # pathway mode: rotated pathway name text above annotation bars
  if (mode == "pathway") {
    pw_rle <- rle(as.character(col_anno$pure_pathway)); pw_ends <- cumsum(pw_rle$lengths)
    pw_starts <- c(1, pw_ends[-length(pw_ends)] + 1); pw_centers <- (pw_starts + pw_ends) / (2 * ncol(beta_mat))
    mat_pos <- p$gtable$layout[grep("matrix", p$gtable$layout$name), ][1, ]; anno_t <- mat_pos$t - 1
    p$gtable <- gtable_add_rows(p$gtable, unit(3.5, "lines"), pos = anno_t - 1)
    txt <- textGrob(label = pw_rle$values, x = pw_centers, y = 0.5, just = c("center","center"),
                    gp = gpar(fontsize = 8, fontface = "bold"), rot = 90)
    p$gtable <- gtable_add_grob(p$gtable, txt, t = anno_t, l = mat_pos$l, r = mat_pos$r, name = "pw_labels", clip = "off")
  }
  
  grid.newpage(); grid.draw(p$gtable)
  return(col_anno)
}


sel_df<-read.xlsx("E:/AID cohort/code/NMF/Spectra/lam0001_dict242_lmem_results_annotated.xlsx")
sel_df

sel_vec <- sel_df[sel_df$celltype=="all", 1]
sel_vec2 <- sel_df[!sel_df$celltype=="all", 1]


navy_blue_yellow <- colorRampPalette(c("#08306B", "#6BAED6", "#FDD835", "#F16913"))(100)
navy_lightblue_orange <- colorRampPalette(c("#0B1F4B", "#4C9ED9", "#BFDDF2", "#FDB863", "#E66101"))(100)
navy_blue_orange <- colorRampPalette(c("#082567", "#2166AC", "#92C5DE", "#FDDDA0", "#D6604D"))(100)
navy_blue_pink_red <- colorRampPalette(c("#08306B", "#6BAED6", "#F4A6B8", "#B2182B"))(100)
navy_blue_red <- colorRampPalette(c("#08306B", "#6BAED6", "#FCA0A0", "#D7191C"))(100)
navy_blue_red2 <- colorRampPalette(c("#053061", "#2166AC", "#F4A6A6", "#B40426"))(100)
navy_blue_red3 <- colorRampPalette(c("#08306B", "#6BAED6", "#FEE0D2", "#CB181D"))(100)

all_groupbycells <- plot_lmem_heatmap(lmem_res, mode = "cells", factor_list = sel_vec, strip_celltype_prefix = TRUE, cellheight = 12,heatmap_cols =navy_blue_red3  )
all_groupbyimmcate<- plot_lmem_heatmap(lmem_res, mode = "immunecate", factor_list = sel_vec, strip_celltype_prefix = TRUE, cellheight = 12,heatmap_cols = navy_blue_red3 )
all_groupbypathway<- plot_lmem_heatmap(lmem_res, mode = "pathway", factor_list = sel_vec, cellheight = 12,heatmap_cols = navy_blue_red3)
cell_groupbycells <- plot_lmem_heatmap(lmem_res, mode = "cells", factor_list = sel_vec2, strip_celltype_prefix = TRUE, cellheight = 12,heatmap_cols = navy_blue_red3)
cell_groupbyimmcate<- plot_lmem_heatmap(lmem_res, mode = "immunecate", factor_list = sel_vec2, strip_celltype_prefix = TRUE, cellheight = 12,heatmap_cols = navy_blue_red3)
cell_groupbypathway<- plot_lmem_heatmap(lmem_res, mode = "pathway", factor_list = sel_vec2, cellheight = 12)


###(2)discreate color values for beta (Finally,use this)
# dis_col1=c("#0D56A0", "#2F80BC", "#6EAFD5", "#B1D2E8","#FCDECC", "#F4AA87", "#D25748", "#B82530")
# dis_col1=c("#464D8A","#6673B5","#8F98C8","#B8B5D6","#D8C8D8","#F6D7D7","#E9A6A6","#C75A62")
# dis_col1=c("#0D56A0", "#2F80BC", "#6EAFD5", "#B1D2E8","#EEDCAC","#D2A75C","#BE812C","#754307")
# dis_col1=c("#464D8A","#6673B5","#8F98C8","#B8B5D6","#FFEAE3","#F4AA87", "#D25748", "#B82530")
# dis_col1=c("#464D8A","#6673B5","#8F98C8","#B8B5D6","#F7E9D9","#ECC79F","#DFA067","#BE812C")
dis_col1=c("#2F4692","#335FAA","#4979BB","#8EADD4","#FFF1E9","#FCD9BB","#E89F60","#BE812C")

plot_lmem_heatmap_binned <- function(lmem_res,
                                     mode = c("cells", "immunecate", "pathway"),
                                     factor_list = NULL,
                                     strip_celltype_prefix = FALSE,
                                     cellwidth = NA,
                                     cellheight = 10,
                                     beta_breaks = c(-Inf, -18, -12, -6, 0, 6, 12, 18, Inf),
                                     beta_cols = dis_col1,
                                     beta_labels = c("<=-18", "-18~-12", "-12~-6", "-6~0",
                                                     "0~6", "6~12", "12~18", ">18"),
                                     immune_rds = "E:\\AID cohort\\code\\NMF\\Spectra\\0830_PieImmuneCate.RDS",
                                     na_text_size = 4,
                                     legend_height = 0.4,
                                     title = NULL) {
  mode <- match.arg(mode)
  if (is.null(title)) title <- sprintf("%s_groupby%s_binned (LMEM beta, disease vs HC)",
                                       ifelse(is.null(factor_list), "whole", "partial"), mode)
  library(pheatmap); library(dplyr); library(tidyr); library(stringr); library(gtable); library(grid)
  
  # matrices
  beta_mat <- lmem_res %>% select(disease, factor, beta) %>%
    pivot_wider(names_from = factor, values_from = beta) %>% column_to_rownames("disease") %>% as.matrix()
  sig_mat <- lmem_res %>% select(disease, factor, fdr_global) %>%
    mutate(sig = ifelse(fdr_global < 0.05, "*", "")) %>% select(-fdr_global) %>%
    pivot_wider(names_from = factor, values_from = sig) %>% column_to_rownames("disease") %>% as.matrix()
  
  # palettes
  ct_order <- c("DC","Mono","pDC","LDG","Baso","B","Plasma",
                "NnCD4T","NonNCD4T","NnCD8T","NonNCD8T","UTC","NK","Platelet","all")
  celltype_colors <- c("DC"="#F0E442","Mono"="#F49D5C","pDC"="#B07A99","LDG"="#C95968",
                       "Baso"="#8C2522","B"="#003D81","Plasma"="#B3D1E7","NnCD4T"="#847AB3","NonNCD4T"="#C8BBD5",
                       "NnCD8T"="#F6E1EE","NonNCD8T"="#0B71AB","UTC"="#86C7B4","NK"="#96AF95","Platelet"="#B3928B","all"="#502680")
  Immune_cateOrder <- c("Cell function and homeostasis","Inflammatory regulation","Adaptive immune regulation",
                        "Innate immune sensing","Amino acid metabolism","Energy metabolism","Lipid metabolism",
                        "Small molecule & nucleotide metabolism","Intracellular signal transduction")
  imm_colors <- setNames(c("#2F91A7","#7DB6C1","#C1DCE2","#8E9FC9","#D5B3CD","#AA94C1","#E9B48F","#E1CF84","#E4935D"), Immune_cateOrder)
  sp_colors <- c("cell-type-specific"="#7777BD", "Non cell-type-specific"="#CF7D83")
  group_colors <- c("polygenic patient"="#B18749", "monogenic patient"="#8C2522")
  disease_colors <- c("KD"="#1976D2","JIA"="#7F7FD5","SLE"="#998EC3","Uncategorized"="#B07AA1",
                      "Protein_homeostasis"="#DC9497","TBK1_IRF3"="#9D7660","Cytoskeleton_and_small_GTPase"="#BAB0AC",
                      "Negative_regulation_of_IFN-I"="#76B7B2","Inborn_errors_of_cell_death"="#F28E2B","Immune_metabolic"="#EDC948",
                      "Inflammasome_IL-1β"="#8A9FD1","Arachidonic acid metabolism"="#61BADA","Endolysosomal_nucleic_acid_sensing"="#D37295",
                      "Ca2+_flux-PLC"="#A0CBE8","Osteoclast function"="#FABFD2","NFKB pathway"="#FFBE7D")
  
  # col_anno
  factor_ct_vec <- setNames(lmem_res$celltype, lmem_res$factor)
  all_factors <- colnames(beta_mat)
  col_anno <- data.frame(factor_id = all_factors, celltype = factor_ct_vec[all_factors],
                         specificity = ifelse(factor_ct_vec[all_factors] == "all", "Non cell-type-specific", "cell-type-specific"),
                         name = sub("^Factor_\\d+[ _]", "", all_factors), row.names = all_factors, stringsAsFactors = FALSE)
  col_anno$pure_pathway <- sub("^[^_]+_", "", col_anno$name)
  ImmuneCate <- readRDS(immune_rds)
  ct_prefix_regex <- paste0("^(", paste(rev(ct_order), collapse="|"), ")_")
  col_anno$immune_category <- ImmuneCate$new_category_en[match(sub(ct_prefix_regex, "", col_anno$name), ImmuneCate$name)]
  
  # partial filter
  if (!is.null(factor_list)) {
    sel <- factor_list[factor_list != ""]
    matched <- col_anno$factor_id %in% sel | col_anno$name %in% sel
    cat(sprintf("[partial] %d / %d factors matched\n", sum(matched), length(sel)))
    col_anno <- col_anno[matched, ]; beta_mat <- beta_mat[, rownames(col_anno)]; sig_mat <- sig_mat[, rownames(col_anno)]
  }
  
  # sort + gaps by mode
  if (mode == "cells") {
    col_anno <- col_anno[order(factor(col_anno$celltype, levels = ct_order)), ]
    ct_counts <- table(factor(col_anno$celltype, levels = ct_order)); ct_counts <- ct_counts[ct_counts > 0]
    gaps_col <- cumsum(ct_counts)[-length(ct_counts)]
  } else if (mode == "immunecate") {
    col_anno <- col_anno[order(factor(col_anno$immune_category, levels = Immune_cateOrder), factor(col_anno$celltype, levels = ct_order)), ]
    imm_counts <- table(factor(col_anno$immune_category, levels = Immune_cateOrder)); imm_counts <- imm_counts[imm_counts > 0]
    gaps_col <- cumsum(imm_counts)[-length(imm_counts)]
  } else {
    pw_order <- names(sort(table(col_anno$pure_pathway), decreasing = TRUE))
    col_anno <- col_anno[order(factor(col_anno$pure_pathway, levels = pw_order), factor(col_anno$celltype, levels = ct_order)), ]
    pw_counts <- table(factor(col_anno$pure_pathway, levels = pw_order)); pw_counts <- pw_counts[pw_counts > 0]
    gaps_col <- cumsum(pw_counts)[-length(pw_counts)]
  }
  beta_mat <- beta_mat[, rownames(col_anno)]; sig_mat <- sig_mat[, rownames(col_anno)]
  colnames(beta_mat) <- col_anno$name; colnames(sig_mat) <- col_anno$name; rownames(col_anno) <- col_anno$name
  col_anno_plot <- col_anno[, c("celltype", "specificity", "immune_category")]
  
  # column labels
  if (mode == "pathway") {
    labels_col <- rep("", nrow(col_anno))
  } else if (strip_celltype_prefix) {
    labels_col <- sub("^[^_]+_", "", colnames(beta_mat))
  } else { labels_col <- colnames(beta_mat) }
  
  # row order + annotation
  condition_order <- c("KD","JIA","SLE","Uncategorized","Protein_homeostasis","TBK1_IRF3",
                       "Cytoskeleton_and_small_GTPase","Negative_regulation_of_IFN-I","Inborn_errors_of_cell_death",
                       "Immune_metabolic","Inflammasome_IL-1β","Arachidonic acid metabolism","Endolysosomal_nucleic_acid_sensing",
                       "Ca2+_flux-PLC","Osteoclast function","NFKB pathway")
  beta_mat <- beta_mat[condition_order, ]; sig_mat <- sig_mat[condition_order, ]
  disease_group <- data.frame(group = c(rep("polygenic patient", 3), rep("monogenic patient", 13)),
                              disease = condition_order, row.names = condition_order)
  
  # annotation colors
  actual_ct <- unique(col_anno$celltype); actual_imm <- unique(col_anno$immune_category); actual_sp <- unique(col_anno$specificity)
  ann_colors <- list(group = group_colors, disease = disease_colors[condition_order],
                     celltype = celltype_colors[names(celltype_colors) %in% actual_ct],
                     specificity = sp_colors[names(sp_colors) %in% actual_sp],
                     immune_category = imm_colors[names(imm_colors) %in% actual_imm])
  
  # cell coloring: FULL 8 levels, unchanged
  bk_full <- beta_breaks; bk_full[bk_full == -Inf] <- -100; bk_full[bk_full == Inf] <- 100
  
  # legend: bins actually covered
  rng <- range(beta_mat, na.rm = TRUE)
  b1 <- findInterval(rng[1], beta_breaks)
  b2 <- findInterval(rng[2], beta_breaks)
  used <- b1:b2
  leg_col <- beta_cols[used]
  leg_lab <- beta_labels[used]
  cat(sprintf("[legend] data [%.2f, %.2f] -> bins %d:%d, %d blocks: %s\n",
              rng[1], rng[2], b1, b2, length(used), paste(leg_lab, collapse=", ")))
  
  n_col <- ncol(beta_mat); fontsize_col <- ifelse(n_col > 150, 5, ifelse(n_col > 100, 6, 7))
  p <- pheatmap(beta_mat,
                color = beta_cols, breaks = bk_full,
                cluster_rows = FALSE, cluster_cols = FALSE,
                gaps_row = 3, gaps_col = gaps_col,
                annotation_col = col_anno_plot, annotation_row = disease_group,
                annotation_colors = ann_colors,
                annotation_names_row = FALSE, show_rownames = FALSE,
                labels_col = labels_col,
                display_numbers = sig_mat, number_color = "black", fontsize_number = 7,
                border_color = NA,
                main = title,
                angle_col = 90, fontsize_row = 10, fontsize_col = fontsize_col,
                legend = TRUE, cellwidth = cellwidth, cellheight = cellheight, silent = TRUE)
  
  # replace legend: compact height, top-aligned, no borders
  leg_idx <- which(p$gtable$layout$name == "legend")
  if (length(leg_idx) > 0) {
    n <- length(leg_col)
    lh <- legend_height
    gl <- gList()
    for (i in seq_len(n)) {
      yc <- 1 - (i - 0.5) / n * lh
      gl <- gList(gl, rectGrob(x = 0.25, y = yc, width = 0.35, height = 0.95/n * lh,
                               just = "center", gp = gpar(fill = leg_col[i], col = NA)))
      gl <- gList(gl, textGrob(leg_lab[i], x = 0.48, y = yc, just = "left", gp = gpar(fontsize = 8)))
    }
    p$gtable$grobs[[leg_idx]] <- gTree(children = gl)
    cat(sprintf("[legend] %d blocks, height=%.2f\n", n, lh))
  }
  
  # shrink annotation legend color blocks
  ann_leg_idx <- grep("annotation_legend", p$gtable$layout$name)
  if (length(ann_leg_idx) > 0) {
    for (cl in unique(p$gtable$layout$l[ann_leg_idx])) p$gtable$widths[cl] <- unit(0.85, "cm")
    for (i in ann_leg_idx) {
      g <- editGrob(p$gtable$grobs[[i]], gPath("rect"), grep = TRUE, global = TRUE,
                    height = unit(0.32, "cm"), width = unit(0.32, "cm"))
      p$gtable$grobs[[i]] <- g
    }
  }
  
  # shrink NA text
  shrink_na <- function(g) {
    if (is.grob(g)) {
      if (inherits(g, "text") && !is.null(g$label) && ("NA" %in% g$label)) g$gp$fontsize <- na_text_size
      if (!is.null(g$children)) g$children <- lapply(g$children, shrink_na)
    }
    return(g)
  }
  p$gtable$grobs <- lapply(p$gtable$grobs, shrink_na)
  
  # pathway mode: rotated text
  if (mode == "pathway") {
    pw_rle <- rle(as.character(col_anno$pure_pathway)); pw_ends <- cumsum(pw_rle$lengths)
    pw_starts <- c(1, pw_ends[-length(pw_ends)] + 1); pw_centers <- (pw_starts + pw_ends) / (2 * ncol(beta_mat))
    mat_pos <- p$gtable$layout[grep("matrix", p$gtable$layout$name), ][1, ]; anno_t <- mat_pos$t - 1
    p$gtable <- gtable_add_rows(p$gtable, unit(3.5, "lines"), pos = anno_t - 1)
    txt <- textGrob(label = pw_rle$values, x = pw_centers, y = 0.5, just = c("center","center"),
                    gp = gpar(fontsize = 8, fontface = "bold"), rot = 90)
    p$gtable <- gtable_add_grob(p$gtable, txt, t = anno_t, l = mat_pos$l, r = mat_pos$r, name = "pw_labels", clip = "off")
  }
  
  grid.newpage(); grid.draw(p$gtable)
  return(col_anno)
}



sel_df  <- read.xlsx("E:/AID cohort/code/NMF/Spectra/lam0001_dict242_lmem_results_annotated.xlsx")
sel_vec  <- sel_df[sel_df$celltype == "all", 1]
sel_df <- read.xlsx("E:/AID cohort/code/NMF/Spectra/selected_leme2.xlsx")
sel_vec2 <- sel_df[!sel_df$celltype == "all", 1]
unique(sel_df$name)

all_groupbyimmcate_binned <- plot_lmem_heatmap_binned(lmem_res, mode = "immunecate", factor_list = sel_vec,  strip_celltype_prefix = TRUE, cellheight = 10,cellwidth = 8.4)

cell_groupbycells_binned   <- plot_lmem_heatmap_binned(lmem_res, mode = "cells",     factor_list = sel_vec2, strip_celltype_prefix = TRUE, cellheight = 10,cellwidth = 4.4)
#cell_groupbyimmcate_binned <- plot_lmem_heatmap_binned(lmem_res, mode = "immunecate", factor_list = sel_vec2, strip_celltype_prefix = TRUE, cellheight = 12)
cell_groupbypathway_binned <- plot_lmem_heatmap_binned(lmem_res, mode = "pathway",   factor_list = sel_vec2, cellheight = 10,cellwidth = 5.63)

########################################################
##ulm_res celltype id dependent of my_dict(celltype)
###### STEP7: merge annotations into lmem_res, output xlsx + txt
library(openxlsx)
anno_cols <- col_anno_cells[, c("factor_id", "name", "specificity", "immune_category")]

lmem_res_out <- lmem_res %>%
  left_join(anno_cols, by = c("factor" = "factor_id")) %>%
  arrange(fdr_global, p_value)

col_order_out <- c("factor", "name", "celltype", "specificity", "immune_category",
                   "disease", "n_hc", "n_dis", "beta", "se", "t_val",
                   "p_value", "fdr_global", "fdr_celltype", "converged", "note")
lmem_res_out <- lmem_res_out[, col_order_out]

cat("merged rows:", nrow(lmem_res_out), "\n")
cat("NA immune_category:", sum(is.na(lmem_res_out$immune_category)), "\n")
cat("significant (fdr<0.05):", sum(lmem_res_out$fdr_global < 0.05, na.rm = TRUE), "\n")

write.xlsx(lmem_res_out, "lam0001_dict242_lmem_results_annotated.xlsx",
           row.names = FALSE, overwrite = TRUE)
write.table(lmem_res_out, "lam0001_dict242_lmem_results_annotated.txt",
            sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
cat("saved: xlsx + txt\n")
