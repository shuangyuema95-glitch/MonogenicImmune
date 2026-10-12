library(Seurat)
library(tidyverse)
library(scales)
pbmc1<-readRDS("E:\\AID cohort\\code\\result\\PBMC1.rds")

gene_list <- list(
  "NFKB pathway" = c("TNFAIP3", "IRAK2", "NOD2", "RELA","IL1R1","TLR1"),
  "Uncategorized" = c("ADA2", "ELF4", "CSF3R", "STAT4"),
  "Osteoclast function" = c("ACP5", "OGFRL1"),
  "Ca2+_flux-PLC" = c("PLCG2", "PLCG1"),
  "Endolysosomal_nucleic_acid_sensing" = c("UNC93B1", "TLR8", "TLR7","PLD4", "TRAF3"),
  "Arachidonic acid metabolism" = c("TBXAS1"),
  "Inflammasome_IL-1β" = c("NLRC4", "NLRP3", "PSTPIP1", "CDC42", "LPIN2"),
  "Immune_metabolic" = c("LACC1", "SLC7A7"),
  "PP" = c("polygenic"),
  "HC" = c("wide-type"),
  "Inborn_errors_of_cell_death" = c("RIPK1", "TNFRSF1A", "OTULIN", "RNF31"),
  "Negative_regulation_of_IFN-I" = c("ISG15", "USP18", "SOCS1"),
  "Cytoskeleton_and_small_GTPase" = c("GNAI2", "KRAS"),
  "TBK1_IRF3" = c("STING", "TREX1", "IFIH1", "COPA"),
  "Protein_homeostasis" = c("UBA1", "PSMD12"))

gene_map <- stack(gene_list)
colnames(gene_map) <- c("gene", "condition0")
pbmc1$condition0 <- gene_map$condition0[match(pbmc1$gene, gene_map$gene)]
pbmc1@meta.data <- pbmc1@meta.data %>% mutate(
  condition = case_when(
    str_detect(samples, "SLE") ~ "SLE",
    str_detect(samples, "JIA") ~ "JIA",
    str_detect(samples, "KD") ~ "KD",
    str_detect(samples, "HC") ~ "HC",
    TRUE ~ condition0
  )
)
table(pbmc1$condition, useNA="ifany")

######1 Level 1 bubble plots  and classfication for porportation patterns
#' Generate Cell-type Proportion Distribution Bubble + Pie Plot
#' @title plot_celltype_distribution
#' @description Visualize cell-type composition across 17 groups (HC + disease conditions).
#' Main panel: x=group, y=cell type (C01-C15), bubble size & color mapped to within-group proportion.
#' Right-side donut pies: for each cell type, distribution of the 17 groups across 4 proportion tiers (<1%, 1-5%, 5-10%, >10%).
#' Cells are classified into 4 bands based on how many groups they are detectable in (prop >= detect_cut):
#'   Universal (n_found >= universal_min), Shared (intermediate), Specialized (n_found <= specialized_max), Rare (n_found == 0).
#' Bands are ordered top-to-bottom: Universal -> Shared -> Specialized -> Rare, within each band sorted by C number ascending.
#' Dashed horizontal lines separate bands; return combined_table for QC.
#' @param seu Seurat object
#' @param group_col Character, metadata column storing condition/group assignment
#' @param level1_col Character, metadata column storing Level1 cell-type annotation, default "Level1"
#' @param group_order Character vector, fixed group order (left-to-right on x axis)
#' @param detect_cut Numeric, proportion threshold for "detectable", default 0.01
#' @param high_cut Numeric, proportion threshold for "high abundance", default 0.10
#' @param universal_min Integer, minimum number of groups to call Universal, default 12 (17*0.7)
#' @param specialized_max Integer, maximum number of groups to call Specialized, default 2
#' @param tier_breaks Numeric vector, tier boundaries for pies, default c(0.01,0.05,0.10)
#' @param tier_labels Character vector, tier labels, default c("<1%","1-5%","5-10%",">10%")
#' @param tier_colors Named character vector, fill colors for pie tiers
#' @param bubble_ramp Character vector, color ramp for continuous bubble fill
#' @param orientation Character, "v" (group on x, cell on y, pies on right) or "h" (cell on x, group on y, pies on top)
#' @returnType list
#' @return list with $plot (ggplot/patchwork) and $combined_table (cell x n_found x n_high x band x 17 groups proportion)
#' @author MSY

library(RColorBrewer)
col2 <- colorRampPalette(c("#BEBEBE",brewer.pal(9,"YlGnBu")))(100)[1:85]


plot_celltype_distribution <- function(seu,
                                       group_col,
                                       level1_col = "Level1",
                                       group_order = c("HC","KD","JIA","SLE","Uncategorized","Protein_homeostasis",
                                                       "TBK1_IRF3","Cytoskeleton_and_small_GTPase","Negative_regulation_of_IFN-I",
                                                       "Inborn_errors_of_cell_death","Immune_metabolic","Inflammasome_IL-1β",
                                                       "Arachidonic acid metabolism","Endolysosomal_nucleic_acid_sensing",
                                                       "Ca2+_flux-PLC","Osteoclast function","NFKB pathway"),
                                       detect_cut = 0.01,
                                       high_cut = 0.10,
                                       universal_min = 12,
                                       specialized_max = 2,
                                       tier_breaks = c(0.01, 0.05, 0.10),
                                       tier_labels = c("<1%","1-5%","5-10%",">10%"),
                                       tier_colors = c("<1%"="#E8F0F8","1-5%"="#9EC9E0","5-10%"="#3D7DB8",">10%"="#0B3D7A"),
                                       bubble_ramp = col2,
                                       orientation = c("v","h")) {
  orientation <- match.arg(orientation)
  lab_map <- c("HC"="HC","KD"="KD","JIA"="JIA","SLE"="SLE",
               "Osteoclast function"="OST","Ca2+_flux-PLC"="CFP",
               "Endolysosomal_nucleic_acid_sensing"="ENS","Arachidonic acid metabolism"="AAM",
               "Inflammasome_IL-1β"="IIL","Immune_metabolic"="IMM",
               "Inborn_errors_of_cell_death"="IECD","Negative_regulation_of_IFN-I"="NRI",
               "Cytoskeleton_and_small_GTPase"="CSG","TBK1_IRF3"="TIR",
               "Protein_homeostasis"="PHS","NFKB pathway"="NFB","Uncategorized"="UNS")
  c_map <- c("Monocyte"="C01","Naive CD4 T"="C02","Non-naive CD8 T"="C03","B cell"="C04",
             "Non-naive CD4 T"="C05","Naive CD8 T"="C06","NK"="C07","LDG"="C08",
             "Platelet"="C09","Erythrocyte"="C10","Plasma"="C11","DC"="C12",
             "UTC"="C13","Basophil"="C14","pDC"="C15")
  meta <- seu@meta.data
  counts <- as.data.frame(table(meta[[group_col]], meta[[level1_col]]))
  names(counts) <- c("group","Level1","n")
  df <- counts %>% filter(group %in% group_order) %>%
    mutate(cell = recode(Level1, !!!c_map)) %>%
    group_by(group) %>% mutate(prop = n/sum(n)) %>% ungroup() %>%
    mutate(group = factor(group, levels=group_order),
           group_short = factor(lab_map[as.character(group)], levels=lab_map[group_order]),
           tier = cut(prop, breaks=c(-Inf,tier_breaks,Inf), labels=tier_labels))
  band_df <- df %>% group_by(cell) %>%
    summarise(n_found=sum(prop>=detect_cut), n_high=sum(prop>=high_cut), .groups="drop") %>%
    mutate(band = case_when(n_found==0                 ~ "Rare",
                            n_found>=universal_min    ~ "Universal",
                            n_found<=specialized_max   ~ "Specialized",
                            TRUE                      ~ "Shared"))
  df <- df %>% left_join(band_df, by="cell")
  ord_df <- band_df %>%
    mutate(band=factor(band,levels=c("Universal","Shared","Specialized","Rare")),
           cnum=as.integer(sub("C","",cell))) %>%
    arrange(band,cnum)
  cat("y-axis top-to-bottom order:\n"); print(ord_df %>% select(cell,band,cnum), n=Inf)
  cell_levels <- rev(ord_df$cell)
  df <- df %>% mutate(cell=factor(cell,levels=cell_levels),
                      band=factor(band,levels=c("Universal","Shared","Specialized","Rare")))
  n_univ <- sum(band_df$band=="Universal"); n_shar <- sum(band_df$band=="Shared")
  n_spec <- sum(band_df$band=="Specialized"); n_rare <- sum(band_df$band=="Rare")
  pie_df <- df %>% count(cell,tier,.drop=FALSE) %>% group_by(cell) %>%
    mutate(freq=n/sum(n)) %>% ungroup()
  prop_mat <- df %>% select(cell, group, prop) %>%
    pivot_wider(names_from=group, values_from=prop) %>%
    arrange(match(cell, cell_levels))
  band_out <- band_df %>% arrange(desc(n_found)) %>% select(cell, n_found, n_high, band)
  combined_table <- band_out %>% left_join(prop_mat, by="cell")
  if (orientation=="v") {
    s1 <- n_rare + 0.5
    s2 <- n_rare + n_spec + 0.5
    s3 <- n_rare + n_spec + n_shar + 0.5
    main <- ggplot(df,aes(x=group_short,y=cell)) +
      geom_point(aes(size=prop,color=prop),alpha=0.85) +
      geom_hline(yintercept=c(s1,s2,s3),linetype="dashed",color="black",linewidth=0.4) +
      annotate("text",x=0.4,y=(n_rare+1)/2,label="Rare",angle=90,hjust=0.5,size=3.5) +
      annotate("text",x=0.4,y=n_rare+(n_spec+1)/2,label="Specialized",angle=90,hjust=0.5,size=3.5) +
      annotate("text",x=0.4,y=n_rare+n_spec+(n_shar+1)/2,label="Shared",angle=90,hjust=0.5,size=3.5) +
      annotate("text",x=0.4,y=n_rare+n_spec+n_shar+(n_univ+1)/2,label="Universal",angle=90,hjust=0.5,size=3.5) +
      scale_size_continuous(range=c(0.6,5),name="Proportion") +
      scale_color_gradientn(colours=bubble_ramp,name="Proportion") +
      theme_classic() +
      theme(panel.border=element_rect(color="black",fill=NA,linewidth=0.6),
            axis.text=element_text(size=9,color="black"),
            axis.text.x=element_text(angle=45,hjust=1,color="black"),
            axis.title=element_blank(),axis.ticks=element_line(color="black"))
    pies <- lapply(ord_df$cell,function(ct)
      ggplot(pie_df%>%filter(cell==ct),aes(x=1,y=freq,fill=tier)) +
        geom_col(width=1,color=NA,linewidth=0) + coord_polar("y") +
        scale_fill_manual(values=tier_colors,drop=FALSE) + theme_void() + theme(legend.position="none"))
    p <- main + wrap_plots(pies,nrow=length(pies)) + plot_layout(widths=c(6,1))
  } else {
    s1 <- n_univ+0.5; s2 <- n_univ+n_shar+0.5; s3 <- n_univ+n_shar+n_spec+0.5
    main <- ggplot(df,aes(x=cell,y=group_short)) +
      geom_point(aes(size=prop,color=prop),alpha=0.85) +
      geom_vline(xintercept=c(s1,s2,s3),linetype="dashed",color="black",linewidth=0.4) +
      annotate("text",y=0.3,x=(n_univ+1)/2,label="Universal",hjust=0.5,vjust=0,size=3.5) +
      annotate("text",y=0.3,x=n_univ+(n_shar+1)/2,label="Shared",hjust=0.5,vjust=0,size=3.5) +
      annotate("text",y=0.3,x=n_univ+n_shar+(n_spec+1)/2,label="Specialized",hjust=0.5,vjust=0,size=3.5) +
      annotate("text",y=0.3,x=n_univ+n_shar+n_spec+(n_rare+1)/2,label="Rare",hjust=0.5,vjust=0,size=3.5) +
      scale_size_continuous(range=c(0.6,5),name="Proportion") +
      scale_color_gradientn(colours=bubble_ramp,name="Proportion") +
      theme_classic() +
      theme(panel.border=element_rect(color="black",fill=NA,linewidth=0.6),
            axis.text=element_text(size=9,color="black"),
            axis.text.x=element_text(angle=45,hjust=1,color="black"),
            axis.title=element_blank(),axis.ticks=element_line(color="black"))
    pies <- lapply(ord_df$cell,function(ct)
      ggplot(pie_df%>%filter(cell==ct),aes(x=1,y=freq,fill=tier)) +
        geom_col(width=1,color=NA,linewidth=0) + coord_polar("y") +
        scale_fill_manual(values=tier_colors,drop=FALSE) + theme_void() + theme(legend.position="none"))
    p <- wrap_plots(pies,ncol=length(pies)) / main + plot_layout(heights=c(1,4))
  }
  return(list(plot=p, combined_table=combined_table))
}

res <- plot_celltype_distribution(pbmc1, group_col="condition", orientation="v",
tier_colors = c("<1%"="#847AB3","1-5%"="#D1929B","5-10%"="#C6D199",">10%"="#0B71AB"))

print(res$plot)
print(res$combined_table, n=Inf, width=Inf)

print(res$plot)
print(res$combined_table, n=Inf, width=Inf)
View(res$combined_table)
library(openxlsx)
write.xlsx(res$combined_table,file="E:\\AID cohort\\code\\CellAnnotation\\Level1_bubble_pie.xlsx")


######1 Level 2 prop stack bar plot (denominator is total count for each level 2 cell)
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
      colour = "white", size = 2.4
    )
  }
  
  return(list(plot = p, prop_table = count_df))
}

gene_list <- list(
  "NFKB pathway" = c("TNFAIP3","IRAK2","NOD2","RELA","TLR1","IL1R1"),
  "Uncategorized" = c("ADA2","ELF4","CSF3R","STAT4"),
  "Osteoclast function" = c("ACP5","OGFRL1"),
  "Ca2+_flux-PLC" = c("PLCG2","PLCG1"),
  "Endolysosomal_nucleic_acid_sensing" = c("UNC93B1","TLR8","TLR7","PLD4","TRAF3"),
  "Arachidonic acid metabolism" = c("TBXAS1"),
  "Inflammasome_IL-1β" = c("NLRC4","NLRP3","PSTPIP1","CDC42","LPIN2"),
  "Immune_metabolic" = c("LACC1","SLC7A7"),
  "Inborn_errors_of_cell_death" = c("RIPK1","TNFRSF1A","OTULIN","RNF31"),
  "Negative_regulation_of_IFN-I" = c("ISG15","USP18","SOCS1"),
  "Cytoskeleton_and_small_GTPase" = c("GNAI2","KRAS"),
  "TBK1_IRF3" = c("STING","TREX1","IFIH1","COPA"),
  "Protein_homeostasis" = c("UBA1","PSMD12"))


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
  show_pct_label = TRUE
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
###Level 2 count as denominator scCODA is table(sample, Level2) raw count as input
