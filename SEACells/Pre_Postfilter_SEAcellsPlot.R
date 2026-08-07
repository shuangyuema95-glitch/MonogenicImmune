library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(cowplot)

csv_dir <- "E:/AID cohort/code/NMF/plot_seacells/"

# lowPurity fill palette
fill_palette <- c("FALSE" = "#003D81", "TRUE" = "#F49D5C")

# Fixed Level1 cell order & color map (space version, match metadata)
fixed_level1_order <- c("DC","Monocyte","pDC","LDG","Basophil","B cell","Plasma","Naive CD4 T","Non-naive CD4 T","Naive CD8 T","Non-naive CD8 T","UTC","NK","Erythrocyte","Platelet")
my_colors <- c(
  "Naive CD4 T"="#847AB3","Non-naive CD4 T"="#C8BBD5",
  "Naive CD8 T"="#F6E1EE","Non-naive CD8 T"="#0B71AB",
  "NK"="#96AF95","B cell"="#003D81","Monocyte"="#F49D5C",
  "LDG"="#C95968","Plasma"="#B3D1E7","Erythrocyte"="#FFEFC1",
  "Platelet"="#B3928B","UTC"="#86C7B4","DC"="#F0E442",
  "pDC"="#B07A99","Basophil"="#8C2222"
)

# Fixed condition x-axis order
condition_order <- c("HC", "KD", "JIA","SLE","Uncategorized","Protein_homeostasis",
                     "TBK1_IRF3","Cytoskeleton_and_small_GTPase","Negative_regulation_of_IFN-I",
                     "Inborn_errors_of_cell_death","Immune_metabolic","Inflammasome_IL-1β",
                     "Arachidonic acid metabolism","Endolysosomal_nucleic_acid_sensing",
                     "Ca2+_flux-PLC","Osteoclast function","NFKB pathway")

# Unified global theme
global_theme <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    text = element_text(color = "black", face = "plain"),
    axis.text = element_text(color = "black", size = 7),
    axis.title = element_text(color = "black", size = 8),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black"),
    plot.title = element_text(color = "black", face = "plain", hjust = 0.5, size = 9),
    legend.text = element_text(color = "black", size = 6),
    legend.title = element_text(color = "black", size = 6),
    legend.key.size = unit(0.3, "cm"),
    panel.border = element_rect(color = "black")
  )

# Helper: fix lowPurity True/False to UPPERCASE
fix_lowpurity_case <- function(df_long){
  df_long$lowPurity <- toupper(df_long$lowPurity)
  return(df_long)
}

# Helper: replace dot back to space for CellType (fix read.csv check.names=T auto replace space to dot)
fix_celltype_dot <- function(df_long){
  df_long$CellType <- gsub("\\.", " ", df_long$CellType)
  return(df_long)
}

###### Pre-filter Plot Functions
plot_ax1 <- function(){
  df <- read.csv(file.path(csv_dir, "ax1_non_trivial_frac.csv"), check.names = FALSE)
  colnames(df) <- c("SEACell_non_trivial_assig", "Fraction_cells")
  ggplot(df, aes(x = factor(SEACell_non_trivial_assig), y = Fraction_cells)) +
    geom_col(fill = "#B3D1E7") +
    labs(title = "Fraction of SEACell_non_trivial_assig",
         x = "SEACell_non_trivial_assig", y = "Fraction cells") +
    global_theme
}

plot_ax2 <- function(){
  df <- read.csv(file.path(csv_dir, "ax2_seacell_cell_count.csv"), check.names = FALSE)
  ggplot(df, aes(x = cell_count)) +
    geom_histogram(bins = 75, fill = "#B3D1E7", color = "white") +
    labs(title = "Number of Cells per SEACell",
         x = "Number of Cells per SEACell", y = "Number of SEACells") +
    global_theme
}

plot_ax3 <- function(){
  df <- read.csv(file.path(csv_dir, "ax3_violin_ct_pre.csv"), check.names = FALSE)
  df$MostAbundantCellType <- factor(df$MostAbundantCellType, levels = fixed_level1_order)
  ggplot(df, aes(x = MostAbundantCellType, y = n_single_cell)) +
    geom_violin(fill = "#B3D1E7", alpha = 0.7) +
    labs(title = "n_single_cell per SEACell (pre‑filter)",
         x = "MostAbundantCellType", y = "n_single_cell") +
    global_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_ax4 <- function(){
  df <- read.csv(file.path(csv_dir, "ax4_ct_mean_prop.csv"), check.names = FALSE)
  df_long <- pivot_longer(df, cols = -MostAbundantCellType, names_to = "CellType", values_to = "Fraction")
  df_long <- fix_celltype_dot(df_long)
  df_long$CellType <- factor(df_long$CellType, levels = fixed_level1_order)
  df_long$MostAbundantCellType <- factor(df_long$MostAbundantCellType, levels = fixed_level1_order)
  ggplot(df_long, aes(x = MostAbundantCellType, y = Fraction, fill = CellType)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = my_colors) +
    labs(title = "Mean cell‑type composition per SEACell (pre‑filter)",
         x = "MostAbundantCellType", y = "Mean fraction") +
    global_theme +
    theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_ax5 <- function(){
  df <- read.csv(file.path(csv_dir, "ax5_ct_type_abs.csv"), check.names = FALSE)
  df$MostAbundantCellType <- factor(df$MostAbundantCellType, levels = fixed_level1_order)
  df_long <- pivot_longer(df, cols = -MostAbundantCellType, names_to = "lowPurity", values_to = "Count")
  df_long <- fix_lowpurity_case(df_long)
  df_long$lowPurity <- factor(df_long$lowPurity, levels = c("FALSE", "TRUE"))
  ggplot(df_long, aes(x = MostAbundantCellType, y = Count, fill = lowPurity)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = fill_palette) +
    labs(x = "MostAbundantCellType", y = "Absolute SEACell number") +
    global_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_ax6 <- function(){
  df <- read.csv(file.path(csv_dir, "ax6_ct_type_norm.csv"), check.names = FALSE)
  df$MostAbundantCellType <- factor(df$MostAbundantCellType, levels = fixed_level1_order)
  df_long <- pivot_longer(df, cols = -MostAbundantCellType, names_to = "lowPurity", values_to = "Freq")
  df_long <- fix_lowpurity_case(df_long)
  df_long$lowPurity <- factor(df_long$lowPurity, levels = c("FALSE", "TRUE"))
  ggplot(df_long, aes(x = MostAbundantCellType, y = Freq, fill = lowPurity)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = fill_palette) +
    labs(x = "MostAbundantCellType", y = "SEACell frequency") +
    global_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_ax7 <- function(){
  df <- read.csv(file.path(csv_dir, "ax7_disease_abs.csv"), check.names = FALSE)
  df_long <- pivot_longer(df, cols = -status, names_to = "lowPurity", values_to = "Count")
  df_long <- fix_lowpurity_case(df_long)
  df_long$lowPurity <- factor(df_long$lowPurity, levels = c("FALSE", "TRUE"))
  ggplot(df_long, aes(x = status, y = Count, fill = lowPurity)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = fill_palette) +
    labs(x = "status", y = "Absolute SEACell number") +
    global_theme
}

plot_ax8 <- function(){
  df <- read.csv(file.path(csv_dir, "ax8_disease_norm.csv"), check.names = FALSE)
  df_long <- pivot_longer(df, cols = -status, names_to = "lowPurity", values_to = "Freq")
  df_long <- fix_lowpurity_case(df_long)
  df_long$lowPurity <- factor(df_long$lowPurity, levels = c("FALSE", "TRUE"))
  ggplot(df_long, aes(x = status, y = Freq, fill = lowPurity)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = fill_palette) +
    labs(x = "status", y = "SEACell frequency") +
    global_theme
}

plot_ax9 <- function(){
  df <- read.csv(file.path(csv_dir, "ax9_condition_abs.csv"), check.names = FALSE)
  df$condition <- factor(df$condition, levels = condition_order)
  df_long <- pivot_longer(df, cols = -condition, names_to = "lowPurity", values_to = "Count")
  df_long <- fix_lowpurity_case(df_long)
  df_long$lowPurity <- factor(df_long$lowPurity, levels = c("FALSE", "TRUE"))
  ggplot(df_long, aes(x = condition, y = Count, fill = lowPurity)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = fill_palette) +
    labs(x = "Condition", y = "Absolute SEACell number") +
    global_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_ax10 <- function(){
  df <- read.csv(file.path(csv_dir, "ax10_condition_norm.csv"), check.names = FALSE)
  df$condition <- factor(df$condition, levels = condition_order)
  df_long <- pivot_longer(df, cols = -condition, names_to = "lowPurity", values_to = "Freq")
  df_long <- fix_lowpurity_case(df_long)
  df_long$lowPurity <- factor(df_long$lowPurity, levels = c("FALSE", "TRUE"))
  ggplot(df_long, aes(x = condition, y = Freq, fill = lowPurity)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = fill_palette) +
    labs(x = "Condition", y = "SEACell frequency") +
    global_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

draw_pre_filter <- function(){
  p1 <- plot_ax1()
  p2 <- plot_ax2()
  p3 <- plot_ax3()
  p4 <- plot_ax4()
  p5 <- plot_ax5()
  p6 <- plot_ax6()
  p7 <- plot_ax7()
  p8 <- plot_ax8()
  p9 <- plot_ax9()
  p10 <- plot_ax10()
  
  full_pre <- (p1 + p2) / (p3 + p4) / (p5 + p6) / (p7 + p8) / (p9 + p10)
  ggsave("seacell_qc_summary_R.pdf", full_pre, width = 18, height = 26, dpi = 150)
  return(list(p1=p1,p2=p2,p3=p3,p4=p4,p5=p5,p6=p6,p7=p7,p8=p8,p9=p9,p10=p10,full_pre=full_pre))
}

######  Post-filter Plot Functions 
plot_post_a <- function(){
  df <- read.csv(file.path(csv_dir, "post_axa_violin_cond.csv"), check.names = FALSE)
  df$condition <- factor(df$condition, levels = condition_order)
  ggplot(df, aes(x = condition, y = n_single_cell)) +
    geom_violin(fill = "#B3D1E7", alpha = 0.7) +
    labs(title = "n_single_cell per SEACell (post‑filter, group by Condition)",
         x = "Condition", y = "n_single_cell") +
    global_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_post_b <- function(){
  df <- read.csv(file.path(csv_dir, "post_axb_violin_ct.csv"), check.names = FALSE)
  df$MostAbundantCellType <- factor(df$MostAbundantCellType, levels = fixed_level1_order)
  ggplot(df, aes(x = MostAbundantCellType, y = n_single_cell)) +
    geom_violin(fill = "#B3D1E7", alpha = 0.7) +
    labs(title = "n_single_cell per SEACell (post‑filter purity >= 0.75)",
         x = "MostAbundantCellType", y = "n_single_cell") +
    global_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_post_c <- function(){
  df <- read.csv(file.path(csv_dir, "post_axc_ct_mean_prop.csv"), check.names = FALSE)
  df_long <- pivot_longer(df, cols = -MostAbundantCellType, names_to = "CellType", values_to = "Fraction")
  df_long <- fix_celltype_dot(df_long)
  df_long$CellType <- factor(df_long$CellType, levels = fixed_level1_order)
  df_long$MostAbundantCellType <- factor(df_long$MostAbundantCellType, levels = fixed_level1_order)
  ggplot(df_long, aes(x = MostAbundantCellType, y = Fraction, fill = CellType)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = my_colors) +
    labs(title = "Mean cell‑type composition per SEACell (post‑filter)",
         x = "MostAbundantCellType", y = "Mean fraction") +
    global_theme +
    theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1))
}

draw_post_filter <- function(){
  pa <- plot_post_a()
  pb <- plot_post_b()
  pc <- plot_post_c()
  
  full_post <- pa / pb / pc
  ggsave("seacell_qc_summary_postfilter_R.pdf", full_post, width = 14, height = 12, dpi = 150)
  return(list(pa=pa,pb=pb,pc=pc,full_post=full_post))
}

# -------------------------- Run Drawing --------------------------
pre_plot_list <- draw_pre_filter()
post_plot_list <- draw_post_filter()

print(pre_plot_list$full_pre)
print(post_plot_list$full_post)

pre_plot_list$p1#no need
pre_plot_list$p2#no need
pre_plot_list$p3#no need
pre_plot_list$p6#no need
pre_plot_list$p7#no need
pre_plot_list$p8#no need
pre_plot_list$p9#no need
post_plot_list$pa#no need
post_plot_list$pb#no need


library(patchwork)
library(cowplot)
p1 <- pre_plot_list$p4 + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_blank(),
        legend.position = "none")

p2 <- post_plot_list$pc + 
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_blank(),
        legend.position = "none")

block_left <- p1 / p2

p3 <- pre_plot_list$p10 + 
  theme(axis.title.x = element_blank(),
        legend.position = "none")

p4 <- pre_plot_list$p6 + 
  theme(axis.title.x = element_blank(),
        legend.position = "none",axis.title.y = element_blank())

block_right <- p3 + p4

final_fig <- block_left | block_right
final_fig <- final_fig + plot_layout(widths = c(1.1, 2.97))
print(final_fig)

ggsave("E:\\AID cohort\\code\\NMF\\plot_seacells\\seacell_custom_layout.pdf", final_fig, width = 9.83, height = 3.75, dpi = 600)
