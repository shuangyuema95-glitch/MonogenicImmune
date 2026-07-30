library(Seurat)
library(data.table)

pbmc1 <- readRDS("/public/home/yuxiaomingroup/msy/work/AID/juzhen/clean/result/50/PBMC1.rds")
gs_df <- fread("0725_FinalList.txt", sep = "\t", header = TRUE)

gene_sig_list <- list()
valid_all_genes <- rownames(pbmc1)

for (i in 1:nrow(gs_df)) {
  gs_name <- gs_df$name[i]
  gene_raw <- strsplit(gs_df$gene[i], ",")[[1]]
  gene_raw <- trimws(gene_raw)
  gene_filter <- intersect(gene_raw, valid_all_genes)
  
  if (length(gene_filter) > 0) {
    gene_sig_list[[gs_name]] <- gene_filter
  } else {
    message(paste0("[SKIP] Pathway ", gs_name, ": zero matched genes"))
  }
}

message(paste0("Total valid pathways to calculate: ", length(gene_sig_list)))

for (sig_nm in names(gene_sig_list)) {
  message(paste0("Start calculating pathway: ", sig_nm))
  pbmc1 <- AddModuleScore(
    object = pbmc1,
    features = list(gene_sig_list[[sig_nm]]),
    nbin = 24,
    ctrl = 100,
    seed = 1,
    assay = "RNA",
    slot = "data",
    name = paste0(sig_nm, "_score")
  )
  message(paste0("Finished pathway: ", sig_nm))
}

saveRDS(pbmc1@meta.data,file="Addmolescore.rds")
meta_df <- pbmc1@meta.data
fwrite(meta_df, file = "PBMC1_AddModuleScore_meta.csv", row.names = TRUE)
message("All calculation finished! Updated Seurat object saved: PBMC1_AddModuleScore_updated.rds")

###Vis
library(data.table)
library(tidyverse)
library(ggplot2)
library(dplyr)
setwd("E:\\AID cohort\\code\\result")
addscore<-read.csv("r13_PBMC1_AddModuleScore_meta.csv")
mapping<-fread("E:\\AID cohort\\code\\ImmuneScoreCal\\scoreMap.txt")
dim(mapping)
colnames(addscore)[grep("_score1",colnames(addscore))]<-mapping$name
colnames(addscore)
unique(mapping$groups)
match(mapping$pathway,colnames(addscore))



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

gene_map <- stack(gene_list)
colnames(gene_map) <- c("gene", "condition0")
addscore$condition0 <- gene_map$condition0[match(addscore$gene, gene_map$gene)]

addscore <- addscore %>% mutate(
  condition = case_when(
    str_detect(samples, "SLE") ~ "SLE",
    str_detect(samples, "JIA") ~ "JIA",
    str_detect(samples, "KD") ~ "KD",
    str_detect(samples, "HC") ~ "HC",
    TRUE ~ condition0
  )
)

unique(addscore$condition)
pathway_names<-mapping$name
agg_long <- addscore %>%
  group_by(Level1, condition) %>%
  summarise(
    across(all_of(pathway_names), ~mean(.x, na.rm = TRUE)),
    .groups = "drop")

agg_long <- agg_long %>%
  pivot_longer(
    cols = all_of(pathway_names),
    names_to = "pathway",
    values_to = "score")

library(tidyverse)

pathway_group_list <- list(
  "Amino acid metabolism" = c(
    "Alanine_Aspartate_Glutamate metabolism",
    "Glycine_Serine_Threonine metabolism",
    "Arginine_and_proline_metabolism",
    "Cysteine_and_methionine_metabolism",
    "Tryptophan_metabolism"
  ),
  "Lipid metabolism" = c(
    "Arachidonic_acid_metabolism",
    "Ether_lipid_metabolism",
    "Fatty_acid_biosynthesis",
    "Fatty_acid_degradation",
    "Glycerolipid_metabolism",
    "Glycerophospholipid_metabolism",
    "Sphingolipid_metabolism"
  ),
  "Glucose metabolism and energy production" = c(
    "Glycolysis",
    "Citrate_cycle_tca_cycle",
    "Oxidative_phosphorylation",
    "Pentose_phosphate_pathway",
    "Pyruvate_metabolism"
  ),
  "Nucleotide metabolism" = c("Purine_metabolism","Pyrimidine_metabolism"),
  "One-carbon and cofactor metabolism" = c("One_carbon_pool_by_folate"),
  "Signaling molecule metabolism" = c("Inositol_phosphate_metabolism","Glutathione_metabolism"),
  "Autophagy and vesicle transport" = c(
    "Autophagy",
    "Mitophagy",
    "CopI_mediated_anterograde_transport",
    "CopII_mediated_vesicle_transport",
    "Endocytosis",
    "Exocytosis"
  ),
  "Protein degradation and homeostasis" = c("Ubiquitin_proteolysis"),
  "Post-translational protein modification" = c("N_glycan_biosynthesis","Unfolded_protein_response"),
  "Oxidative stress and cellular stress" = c("Oxidative_stress","Hypoxia"),
  "Cell death" = c(
    "Pro_apoptosis",
    "Anti_apoptosis",
    "Ferroptosis",
    "Necroptosis",
    "Netosis",
    "Pyroptosis",
    "Defective_RIPK1_mediated_regulated_necrosis",
    "TRAIL_signaling",
    "Inflammasomes"
  ),
  "Cellular senescence and secretory phenotype" = c("Sasp"),
  "Cell cycle and proliferation" = c("E2F_targets"),
  "Innate immune recognition and signaling" = c(
    "CLR_signaling",
    "NLR_signaling",
    "RLR_signaling",
    "TLR4_signaling_and_tolerance",
    "Toll_like_receptor_9_TLR9_cascade",
    "Toll_like_receptor_TLR1_TLR2_cascade",
    "MYD88_distinct_inputoutput_pathway"
  ),
  "Intracellular signal transduction" = c(
    "Ca2_PLC_signaling",
    "ERK_MAPK_targets",
    "IRAK2_mediated_activation_of_TAK1_complex",
    "IRAK4_deficiency_TLR2_4",
    "MAPK_signaling_pathway",
    "MTORC1_signaling",
    "NFAT_signaling",
    "NFkB_signaling",
    "PI3K_AKT_MTOR_signaling",
    "TNFA_signaling_via_NFKB",
    "FLT3_signaling",
    "Wnt_beta_catenin_signaling"
  ),
  "Interferon signaling" = c(
    "Interferon_stimulated_genes",
    "IRF3_mediated_induction_of_type_I_IFN",
    "Type_I_interferon_signaling",
    "Type_II_interferon_signaling",
    "ZBP1_DAI_mediated_induction_of_type_I_IFNs",
    "CGAS_STING_signaling"
  ),
  "Cytokine and inflammatory signaling" = c(
    "Chemokine",
    "Chemokine_signaling",
    "Chemokine_receptors",
    "Cytokine_receptors_antiinflammatory",
    "Cytokine_receptors_proinflammatory",
    "IL1_signaling_pathway",
    "IL10_antiinflammatory_signaling_pathway",
    "IL12_IL23_signaling",
    "IL17_signaling",
    "IL18_signaling_pathway",
    "IL2_STAT5_signaling",
    "IL4_IL13_response",
    "IL6_JAK_STAT3_signaling",
    "TNF_family_signaling",
    "TGFb_signaling"
  ),
  "Complement and phagocytosis" = c("Complement_activation","Phagocytosis","Efferocytosis"),
  "Innate immune effector and antimicrobial immunity" = c("Respiratory_burst","Neutrophil_degranulation","Mast_granule_exocytosis"),
  "Myeloid cell development, polarization and function" = c(
    "M1_macrophage",
    "M2_macrophage",
    "Mac_angiogenic_effectors",
    "Mac_CSF1_response",
    "Dc_antigen_crosspresentation",
    "Dc_lps_response",
    "Dendritic_cell_migration",
    "Pdc_cpg_TLR9_response"
  ),
  "B cell development and function" = c("B_cell_receptor_signaling"),
  "Mucosal and humoral immunity" = c("Intestinal_IgA_network","VDJ_recombination"),
  "T cell activation, differentiation and immune checkpoint" = c(
    "CTLA4_inhibitory_signaling",
    "T_cell_receptor_signaling",
    "Th1_differentiation",
    "Th17_differentiation",
    "Th2_differentiation",
    "Treg_differentiation",
    "Treg_FOXP3_stabilization",
    "CD4_Th17",
    "CD8_terminal_exhaustion",
    "Naïve_function",
    "Effector_activation",
    "Anergy"
  ),
  "Lymphocyte development, activation and effector function" = c("Lymphocyte_activation","Cytotoxicity"),
  "Innate lymphoid cell function" = c("NK_cell_cytotoxic"),
  "Immune memory and immunological synapse" = c("Immunological_memory","Immunological_synapse"),
  "Cell adhesion, migration and cytoskeleton regulation" = c(
    "Cell_adhesion",
    "Leukocyte_transendothelial_migration",
    "Lymphocyte_trafficking",
    "Inactivation_of_CDC42_and_RAC1",
    "RAC1_gtpase_cycle",
    "RAC2_gtpase_cycle",
    "RAC3_gtpase_cycle",
    "RHOA_gtpase_cycle",
    "RHOB_gtpase_cycle"
  ),
  "Antigen processing and presentation" = c("MHC_class_I_antigen_presentation","MHC_class_II_antigen_presentation"),
  "Coagulation and hemostasis" = c("Hemostasis"),
  "Skeletal system signaling" = c("Osteoblast_signaling","Osteoclast_signaling")
)

intersect(as.character(unlist(pathway_group_list)),mapping$name)


mapping <- stack(pathway_group_list) %>%
  as_tibble() %>%
  rename(pathway = values, groups = ind)



group_order <- c(
  "Amino acid metabolism",
  "Lipid metabolism",
  "Glucose metabolism and energy production",
  "Nucleotide metabolism",
  "One-carbon and cofactor metabolism",
  "Signaling molecule metabolism",
  "Autophagy and vesicle transport",
  "Protein degradation and homeostasis",
  "Post-translational protein modification",
  "Oxidative stress and cellular stress",
  "Cell death",
  "Cellular senescence and secretory phenotype",
  "Cell cycle and proliferation",
  "Innate immune recognition and signaling",
  "Intracellular signal transduction",
  "Interferon signaling",
  "Cytokine and inflammatory signaling",
  "Complement and phagocytosis",
  "Innate immune effector and antimicrobial immunity",
  "Myeloid cell development, polarization and function",
  "B cell development and function",
  "Mucosal and humoral immunity",
  "T cell activation, differentiation and immune checkpoint",
  "Lymphocyte development, activation and effector function",
  "Innate lymphoid cell function",
  "Immune memory and immunological synapse",
  "Cell adhesion, migration and cytoskeleton regulation",
  "Antigen processing and presentation",
  "Coagulation and hemostasis",
  "Skeletal system signaling"
)

mapping_sorted <- mapping %>%
  mutate(groups = factor(groups, levels = group_order)) %>%
  arrange(groups)

fixed_pathway_order <- mapping_sorted$pathway

length(unlist(pathway_group_list))
setdiff(unique(agg_long$pathway), mapping$pathway)
setdiff(mapping$pathway, unique(agg_long$pathway))

agg_long_anno <- agg_long %>%
  left_join(mapping, by = "pathway")


###(1) plot heatmaps for total pathways
library(tidyverse)
library(pheatmap)
library(ggplotify)
library(patchwork)

plotPathwayHeatmap_A <- function(agg_long, pathway_vec, level1_order, condition_order,
                                 out_pdf, pdf_w=16, pdf_h=14){
  heat_colors <- colorRampPalette(c("#2166ac","white","#b2182b"))(100)
  gg_list <- list()
  for(pw in pathway_vec){
    message("Processing pathway: ", pw)
    pw_dat <- agg_long %>% filter(pathway == pw) %>% select(Level1, condition, score)
    if(nrow(pw_dat)==0){message(pw," skip: no data");next}
    mat_raw <- pivot_wider(pw_dat, names_from=condition, values_from=score, values_fill=NA)
    mat <- mat_raw %>% column_to_rownames("Level1") %>% as.data.frame() %>% mutate(across(everything(),as.numeric))
    mat[is.na(mat)] <- 0
    valid_row <- intersect(level1_order, rownames(mat))
    valid_col <- intersect(condition_order, colnames(mat))
    mat <- mat[valid_row, valid_col, drop=F]
    if(nrow(mat)<2 || ncol(mat)<2 || all(mat==0)){message(pw," skip: insufficient data");next}
    mat_z <- t(scale(t(mat)))
    p <- pheatmap(mat_z, cluster_rows=F, cluster_cols=F, show_rownames=T, show_colnames=T,
                  scale="none", color=heat_colors, main=pw, silent=T)
    gg_obj <- as.ggplot(p)
    gg_list[[length(gg_list)+1]] <- gg_obj
  }
  while(dev.cur() > 1){dev.off()}
  pdf(out_pdf, width=pdf_w, height=pdf_h)
  for(i in seq(1,length(gg_list),by=4)){
    idx <- i:min(i+3,length(gg_list))
    wrap_plots(gg_list[idx],nrow=2,ncol=2) |> print()
  }
  dev.off()
  message("Finished: ",out_pdf)
  return(gg_list)
}

plotPathwayHeatmap_B <- function(agg_long, pathway_vec, level1_order, condition_order,
                                 out_pdf, pdf_w=16, pdf_h=14){
  heat_colors <- colorRampPalette(c("#2166ac","white","#b2182b"))(100)
  gg_list <- list()
  for(pw in pathway_vec){
    message("Processing pathway: ", pw)
    pw_dat <- agg_long %>% filter(pathway == pw) %>% select(Level1, condition, score)
    if(nrow(pw_dat)==0){message(pw," skip: no data");next}
    mat_raw <- pivot_wider(pw_dat, names_from=Level1, values_from=score, values_fill=NA)
    mat <- mat_raw %>% column_to_rownames("condition") %>% as.data.frame() %>% mutate(across(everything(),as.numeric))
    mat[is.na(mat)] <- 0
    valid_row <- intersect(condition_order, rownames(mat))
    valid_col <- intersect(level1_order, colnames(mat))
    mat <- mat[valid_row, valid_col, drop=F]
    if(nrow(mat)<2 || ncol(mat)<2 || all(mat==0)){message(pw," skip: insufficient data");next}
    mat_z <- t(scale(t(mat)))
    p <- pheatmap(mat_z, cluster_rows=F, cluster_cols=F, show_rownames=T, show_colnames=T,
                  scale="none", color=heat_colors, main=pw, silent=T)
    gg_obj <- as.ggplot(p)
    gg_list[[length(gg_list)+1]] <- gg_obj
  }
  while(dev.cur() > 1){dev.off()}
  pdf(out_pdf, width=pdf_w, height=pdf_h)
  for(i in seq(1,length(gg_list),by=4)){
    idx <- i:min(i+3,length(gg_list))
    wrap_plots(gg_list[idx],nrow=2,ncol=2) |> print()
  }
  dev.off()
  message("Finished: ",out_pdf)
  return(gg_list)
}


level1_order <- c("DC","Monocyte","pDC","LDG","Basophil","B cell",
                  "Plasma","Naive CD4 T","Non-naive CD4 T","Naive CD8 T","Non-naive CD8 T",
                  "UTC","NK","Erythrocyte","Platelet")

condition_order <- c("HC", "KD", "JIA","SLE","Uncategorized","Protein_homeostasis",
                     "TBK1_IRF3","Cytoskeleton_and_small_GTPase","Negative_regulation_of_IFN-I","Inborn_errors_of_cell_death",
                     "Immune_metabolic","Inflammasome_IL-1β","Arachidonic acid metabolism",
                     "Endolysosomal_nucleic_acid_sensing",
                     "Ca2+_flux-PLC","Osteoclast function","NFKB pathway")

fixed_pathway_order <- mapping$pathway

gglist_A <- plotPathwayHeatmap_A(agg_long = agg_long_anno,
                                 pathway_vec = fixed_pathway_order,
                                 level1_order = level1_order,
                                 condition_order = condition_order,
                                 out_pdf = "Pathway_CellY_CondX.pdf")

gglist_B <- plotPathwayHeatmap_B(agg_long = agg_long_anno,
                                 pathway_vec = fixed_pathway_order,
                                 level1_order = level1_order,
                                 condition_order = condition_order,
                                 out_pdf = "Pathway_CondY_CellX.pdf")

#pattern identification
library(tidyverse)
library(openxlsx)
library(tidyverse)

detectPathwayPattern <- function(agg_long){
  group_def <- list(
    HC = "HC",
    PP = c("SLE","JIA","KD")
  )
  group_def$MP <- setdiff(unique(agg_long$condition), c(group_def$HC, group_def$PP))
  
  agg_condition <- agg_long %>%
    group_by(pathway, condition) %>%
    summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
    filter(!is.na(mean_score))
  
  res_list <- list()
  pw_vec <- unique(agg_condition$pathway)
  
  for(pw in pw_vec){
    df_pw <- agg_condition %>% filter(pathway == pw)
    rank_df <- df_pw %>% arrange(desc(mean_score))
    rank_str <- paste(rank_df$condition, collapse = ">")
    hc_rank_pos <- which(rank_df$condition == "HC")
    
    # extract group value vector
    get_vec <- function(gname){
      g_cond <- group_def[[gname]]
      df_pw %>% filter(condition %in% g_cond) %>% pull(mean_score)
    }
    v_hc <- get_vec("HC")
    v_pp <- get_vec("PP")
    v_mp <- get_vec("MP")
    
    pattern_tag <- NA
    if(length(v_hc)>0 && length(v_pp)>0 && length(v_mp)>0){
      hc_val <- v_hc
      # All MP > all PP & all PP > HC
      if(all(v_mp > v_pp) && all(v_pp > hc_val)){
        pattern_tag <- "MP>PP>HC"
      }else if(all(v_mp > hc_val) && all(hc_val > v_pp)){
        pattern_tag <- "MP>HC>PP"
      }else if(all(v_pp > v_mp) && all(v_mp > hc_val)){
        pattern_tag <- "PP>MP>HC"
      }else if(all(v_pp > hc_val) && all(hc_val > v_mp)){
        pattern_tag <- "PP>HC>MP"
      }else if(all(hc_val > v_mp) && all(v_mp > v_pp)){
        pattern_tag <- "HC>MP>PP"
      }else if(all(hc_val > v_pp) && all(v_pp > v_mp)){
        pattern_tag <- "HC>PP>MP"
      }
    }
    
    if(is.na(pattern_tag)){
      fn_extract <- function(gname, topn, botn){
        g_cond <- group_def[[gname]]
        sub <- df_pw %>% filter(condition %in% g_cond) %>% arrange(desc(mean_score))
        maxn <- head(sub$condition, topn)
        minn <- tail(sub$condition, botn)
        list(max=maxn, min=minn)
      }
      pp_info <- fn_extract("PP", topn=1, botn=1)
      mp_info <- fn_extract("MP", topn=2, botn=1)
      
      str_pp <- paste0("PP(Max:",paste(pp_info$max,collapse=","),";Min:",paste(pp_info$min,collapse=","),")")
      str_mp <- paste0("MP(Max:",paste(mp_info$max,collapse=","),";Min:",paste(mp_info$min,collapse=","),")")
      pattern_tag <- paste0(str_pp,",",str_mp,",HC_rank=",hc_rank_pos)
    }
    
    res_list[[length(res_list)+1]] <- tibble(
      pathway = pw,
      condition_rank = rank_str,
      pattern = pattern_tag
    )
  }
  bind_rows(res_list)
}
pattern_df <- as.data.frame(detectPathwayPattern(agg_long_anno))
pattern_df$groups<-as.character(as.matrix(mapping[match(pattern_df$pathway,mapping$pathway),'groups']))
write.xlsx(pattern_df,file="E:\\AID cohort\\code\\ImmuneScoreCal\\AddM_AllpathwayPattern.xlsx")

###(2) For specific pathway
pathway_cell_filter <- list(
  "NK_cell_cytotoxic" = c("NK","Non-naive CD8 T"),
  "MHC_class_II_antigen_presentation" = c("DC","pDC","Naive CD4 T","Non-naive CD4 T"),
  "Th17_differentiation" = c("Non-naive CD4 T","Naive CD4 T"),
  "Th1_differentiation" = c("Non-naive CD4 T","Naive CD4 T"),
  "Anergy" = c("Naive CD4 T","Naive CD8 T","Non-naive CD4 T","Non-naive CD8 T"),
  
  "TLR4_signaling_and_tolerance" = c("DC","Monocyte","pDC","LDG"),
  "Toll_like_receptor_9_TLR9_cascade" = c("pDC","B cell","Plasma"),
  "Toll_like_receptor_TLR1_TLR2_cascade" = c("DC","Monocyte","pDC"),
  "NLR_signaling" = c("DC","Monocyte","LDG"),
  
  "Pyroptosis" = c("Plasma","LDG","Basophil"),
  "Necroptosis" = c("DC","Monocyte","pDC","LDG","Naive CD4 T","Naive CD8 T","Non-naive CD4 T","Non-naive CD8 T"),
  "E2F_targets" = c("DC","Monocyte"),
  
  "MAPK_signaling_pathway" = c("DC","Monocyte","pDC","LDG","Basophil","B cell","Plasma","Naive CD4 T","Non-naive CD4 T","Naive CD8 T","Non-naive CD8 T","UTC","NK"),
  "TNFA_signaling_via_NFKB" = c("DC","Monocyte","pDC","LDG","Basophil","B cell","Plasma","Naive CD4 T","Non-naive CD4 T","Naive CD8 T","Non-naive CD8 T","UTC","NK"),
  "Interferon_stimulated_genes" = c("DC","Monocyte","pDC","LDG","Basophil","B cell","Plasma","Naive CD4 T","Non-naive CD4 T","Naive CD8 T","Non-naive CD8 T","UTC","NK"),
  "Ca2_PLC_signaling" = c("Plasma","B cell","Monocyte"),
  "Glycolysis" = c("DC","Monocyte","pDC","LDG","Basophil","B cell","Plasma","Naive CD4 T","Non-naive CD4 T","Naive CD8 T","Non-naive CD8 T","UTC","NK")
  )
idx<-match(names(pathway_cell_filter),fixed_pathway_order)

library(tidyverse)
library(pheatmap)
library(ggplotify)
library(patchwork)

plotPathwayHeatmap_A <- function(agg_long,pathway_vec,level1_order,condition_order,heat_colors=colorRampPalette(c("#2166ac","white","#b2182b"))(100),nrow=2,ncol=2){
  gg_list <- list()
  for(pw in pathway_vec){
    message("Processing: ",pw)
    mat <- agg_long %>% filter(pathway==pw) %>% select(Level1,condition,score) %>% pivot_wider(names_from=condition,values_from=score) %>% column_to_rownames("Level1") %>% as.matrix()
    mat[is.na(mat)] <- 0
    mat <- mat[intersect(level1_order,rownames(mat)),intersect(condition_order,colnames(mat)),drop=FALSE]
    if(nrow(mat)<2||ncol(mat)<2||all(mat==0)){message(pw," skip");next}
    mat_z <- t(scale(t(mat)))
    mat_z[is.nan(mat_z)] <- 0
    p <- pheatmap(mat_z,cluster_rows=FALSE,cluster_cols=FALSE,show_rownames=TRUE,show_colnames=FALSE,fontsize_row=7,scale="none",color=heat_colors,legend=FALSE,border_color=NA,main=pw,fontsize=8,fontsize_number=8,silent=TRUE)
    gg_list[[pw]] <- as.ggplot(p)
  }
  combined <- wrap_plots(gg_list,nrow=nrow,ncol=ncol)
  return(list(individual=gg_list,combined=combined))
}

plotPathwayHeatmap_B <- function(agg_long,pathway_vec,level1_order,condition_order,heat_colors=colorRampPalette(c("#2166ac","white","#b2182b"))(100),nrow=2,ncol=2){
  gg_list <- list()
  for(pw in pathway_vec){
    message("Processing: ",pw)
    mat <- agg_long %>% filter(pathway==pw) %>% select(Level1,condition,score) %>% pivot_wider(names_from=Level1,values_from=score) %>% column_to_rownames("condition") %>% as.matrix()
    mat[is.na(mat)] <- 0
    mat <- mat[intersect(condition_order,rownames(mat)),intersect(level1_order,colnames(mat)),drop=FALSE]
    if(nrow(mat)<2||ncol(mat)<2||all(mat==0)){message(pw," skip");next}
    mat_z <- t(scale(t(mat)))
    mat_z[is.nan(mat_z)] <- 0
    p <- pheatmap(mat_z,cluster_rows=FALSE,cluster_cols=FALSE,show_rownames=FALSE,show_colnames=TRUE,fontsize_col=7,scale="none",color=heat_colors,legend=FALSE,border_color=NA,main=pw,fontsize=8,fontsize_number=8,silent=TRUE)
    gg_list[[pw]] <- as.ggplot(p)
  }
  combined <- wrap_plots(gg_list,nrow=nrow,ncol=ncol)
  return(list(individual=gg_list,combined=combined))
}

# Call
library(RColorBrewer)
library(colorspace)
target_pathway <- names(pathway_cell_filter)
target_pathway <- factor(target_pathway,levels=names(pathway_cell_filter)) %>% as.character()
heat_colors <- colorRampPalette(c("#2166ac","white","#b2182b"))(100)
heat_colors <- colorRampPalette(
  rev(brewer.pal(11, "RdYlBu"))
)(100)
heat_colors <- colorRampPalette(c("#2c0042","white","#e66101"))(100)
res_A <- plotPathwayHeatmap_A(agg_long=agg_long_anno,pathway_vec=target_pathway,level1_order=level1_order,condition_order=condition_order,heat_colors=heat_colors,nrow=3,ncol=6)
res_A$combined

res_B <- plotPathwayHeatmap_B(agg_long=agg_long_anno,pathway_vec=target_pathway,level1_order=level1_order,condition_order=condition_order,heat_colors=heat_colors,nrow=3,ncol=6)
res_B$combined

###(3) plot score for Monogenic and HC, Monogenic status stratified by mutation function
library(tidyverse)
library(patchwork)
library(cowplot)

heat_col <- colorRampPalette(c("#FFFFFF","#E8F0FC","#C2D4F2","#83A9E2","#4082D8","#2A69BD","#3968B9","#4A67B8"))(100)
mech_fill_color <- c("gain-of-function"="#8C2522","deficiency"="#003D81","haploinsufficiency"="#C8BBD5","hypomorphic"="#96AF95")
mech_full_order <- c("gain-of-function","deficiency","haploinsufficiency","hypomorphic")
mech_3keep <- c("gain-of-function","deficiency","haploinsufficiency")
mech_short_map <- c("gain-of-function"="GOF","deficiency"="DE","haploinsufficiency"="HAPLO","hypomorphic"="HYPOR")
heat_y_level <- c("GOF","DE","HAPLO")

select_pathways <- c("NK_cell_cytotoxic","MHC_class_II_antigen_presentation","Th17_differentiation","Th1_differentiation","Anergy","TLR4_signaling_and_tolerance","Toll_like_receptor_9_TLR9_cascade","Toll_like_receptor_TLR1_TLR2_cascade","NLR_signaling","Pyroptosis","Necroptosis","E2F_targets","MAPK_signaling_pathway","TNFA_signaling_via_NFKB","Interferon_stimulated_genes","Ca2_PLC_signaling","Glycolysis")

addscore_filter <- addscore %>% filter(!Level1 %in% c("Erythrocyte","Platelet"), gene != "polygenic")
gene_label_df_raw <- addscore_filter %>% select(gene,mechanism) %>% distinct()
wt_gene <- "wide-type"
df_wt <- gene_label_df_raw %>% filter(gene == wt_gene)
df_other <- gene_label_df_raw %>% filter(gene != wt_gene) %>% mutate(mechanism=factor(mechanism,levels=mech_full_order)) %>% arrange(mechanism)
gene_label_df <- bind_rows(df_wt, df_other) %>% mutate(comb_label=gene)
gene_fixed_order <- gene_label_df$comb_label
valid_pw <- intersect(select_pathways, mapping$pathway)
cell_type_vec <- unique(addscore_filter$Level1)

stat_res_all <- tibble()
plot_list <- list()
fig_batch <- list()

test_cont_table <- function(mat){
  chitest <- chisq.test(mat)
  min_exp <- min(chitest$expected)
  if(min_exp >=5){
    pval <- chitest$p.value
    dfv <- chitest$parameter
  }else{
    ft <- fisher.test(mat, simulate.p.value=TRUE, B=20000)
    pval <- ft$p.value
    dfv <- NA
  }
  tibble(p=pval, df=dfv)
}

plotSinglePathSingleCell <- function(pw, ct, df, label_order, label_meta){
  plot_df <- df %>% filter(Level1==ct) %>% select(gene,mechanism,all_of(pw)) %>% rename(score=all_of(pw)) %>% left_join(label_meta %>% select(gene,mechanism,comb_label), by=c("gene","mechanism")) %>% mutate(comb_label=factor(comb_label,levels=label_order)) %>% filter(!is.na(score))
  if(nrow(plot_df)<3) return(list(fig=NULL,stat_row=NULL))
  if(length(unique(plot_df$comb_label))<2) return(list(fig=NULL,stat_row=NULL))
  
  wt_mean_val <- plot_df %>% filter(gene=="wide-type") %>% pull(score) %>% mean(na.rm=TRUE)
  if(is.na(wt_mean_val)) return(list(fig=NULL,stat_row=NULL))
  kw <- kruskal.test(score ~ comb_label, data=plot_df)
  kw_p <- kw$p.value
  
  cell_cont_raw <- plot_df %>% filter(mechanism %in% mech_3keep) %>% mutate(direction=ifelse(score>wt_mean_val,"UP","DN")) %>% count(mechanism,direction) %>% pivot_wider(names_from=direction,values_from=n,values_fill=0)
  if(nrow(cell_cont_raw)<3) return(list(fig=NULL,stat_row=NULL))
  mat_cell <- column_to_rownames(cell_cont_raw,"mechanism") %>% as.matrix()
  res_cell <- test_cont_table(mat_cell)
  
  gene_group_stat <- plot_df %>% group_by(gene,mechanism) %>% summarise(mean_score=mean(score,na.rm=TRUE),.groups="drop") %>% filter(mechanism %in% mech_3keep) %>% mutate(direction=ifelse(mean_score>wt_mean_val,"UP","DN"))
  gene_cont_raw <- gene_group_stat %>% count(mechanism,direction) %>% pivot_wider(names_from=direction,values_from=n,values_fill=0)
  if(nrow(gene_cont_raw)<3) return(list(fig=NULL,stat_row=NULL))
  mat_gene <- column_to_rownames(gene_cont_raw,"mechanism") %>% as.matrix()
  res_gene <- test_cont_table(mat_gene)
  
  # Count extraction
  get_count <- function(df,mech,dir){
    if(mech %in% rownames(df) && dir %in% colnames(df)) return(df[mech,dir])
    return(0)
  }
  gof_dn <- get_count(mat_gene,"gain-of-function","DN")
  gof_up <- get_count(mat_gene,"gain-of-function","UP")
  def_dn <- get_count(mat_gene,"deficiency","DN")
  def_up <- get_count(mat_gene,"deficiency","UP")
  haplo_dn <- get_count(mat_gene,"haploinsufficiency","DN")
  haplo_up <- get_count(mat_gene,"haploinsufficiency","UP")
  
  # Extract DN gene names separated by ;
  get_gene_list <- function(df,mech_name){
    vec <- df %>% filter(mechanism == mech_name, direction == "DN") %>% pull(gene)
    paste0(vec, collapse = ";")
  }
  gof_dn_genes <- get_gene_list(gene_group_stat,"gain-of-function")
  def_dn_genes <- get_gene_list(gene_group_stat,"deficiency")
  haplo_dn_genes <- get_gene_list(gene_group_stat,"haploinsufficiency")
  
  stat_row <- tibble(
    pathway=pw,
    celltype=ct,
    kruskal_p=kw_p,
    cell_level_p=res_cell$p,
    gene_group_level_p=res_gene$p,
    df_cell=res_cell$df,
    df_genegroup=res_gene$df,
    GOF_DN = gof_dn,
    GOF_UP = gof_up,
    DEF_DN = def_dn,
    DEF_UP = def_up,
    HAPLO_DN = haplo_dn,
    HAPLO_UP = haplo_up,
    GOF_DN_Genes = gof_dn_genes,
    DEF_DN_Genes = def_dn_genes,
    HAPLO_DN_Genes = haplo_dn_genes
  )
  
  heat_df_cell <- cell_cont_raw %>% mutate(mechanism_short=mech_short_map[mechanism]) %>% pivot_longer(-c(mechanism,mechanism_short),names_to="dir",values_to="count") %>% mutate(mechanism_short=factor(mechanism_short,levels=heat_y_level))
  p_heat1 <- ggplot(heat_df_cell,aes(x=dir,y=mechanism_short,fill=count))+geom_tile()+geom_text(aes(label=count),size=2.7)+scale_fill_gradientn(colors=heat_col)+scale_y_discrete(limits=heat_y_level)+labs(title="Cell count")+theme_void()+theme(plot.title=element_text(size=9,hjust=0.5),legend.position="none",axis.text.y=element_text(size=7),axis.text.x=element_text(size=7))
  
  heat_df_gene <- gene_cont_raw %>% mutate(mechanism_short=mech_short_map[mechanism]) %>% pivot_longer(-c(mechanism,mechanism_short),names_to="dir",values_to="count") %>% mutate(mechanism_short=factor(mechanism_short,levels=heat_y_level))
  p_heat2 <- ggplot(heat_df_gene,aes(x=dir,y=mechanism_short,fill=count))+geom_tile()+geom_text(aes(label=count),size=2.7)+scale_fill_gradientn(colors=heat_col)+scale_y_discrete(limits=heat_y_level)+labs(title="Gene group count")+theme_void()+theme(plot.title=element_text(size=9,hjust=0.5),legend.position="none",axis.text.y=element_text(size=7),axis.text.x=element_text(size=7))
  
  heat_combine <- p_heat1 + p_heat2 + plot_layout(ncol=2)
  
  p_box <- ggplot(plot_df,aes(x=comb_label,y=score,fill=mechanism))+
    geom_boxplot(color="black",width=0.6,outlier.shape=NA,linewidth=0.3)+
    geom_hline(yintercept=wt_mean_val,linetype=2,color="gray60",linewidth=0.7)+
    scale_fill_manual(values=mech_fill_color)+
    labs(title=paste0(pw," | ",ct),y="Module Score",x=NULL)+
    theme_classic()+
    theme(axis.text.x=element_text(angle=90,hjust=0.5,colour="black",size=8),axis.text.y=element_text(colour="black",size=8),axis.title=element_text(colour="black"),plot.title=element_text(colour="black",size=10),legend.position="none",axis.line=element_line(color="black"),axis.ticks=element_line(color="black"))
  
  combined <- ggdraw(p_box)+draw_plot(heat_combine,x=0.62,y=0.65,width=0.35,height=0.32)
  list(fig=combined,stat_row=stat_row)
}

while(dev.cur()>1){dev.off()}

for(pw in valid_pw){
  for(ct in cell_type_vec){
    res <- tryCatch(plotSinglePathSingleCell(pw=pw,ct=ct,df=addscore_filter,label_order=gene_fixed_order,label_meta=gene_label_df),error=function(e) list(fig=NULL,stat_row=NULL))
    fig <- res$fig
    statrow <- res$stat_row
    if(!is.null(fig)&&!is.null(statrow)){
      key_name <- paste0(pw,"___",str_replace_all(ct," ","_"))
      plot_list[[key_name]] <- fig
      stat_res_all <- bind_rows(stat_res_all,statrow)
      fig_batch[[length(fig_batch)+1]] <- fig
      if(length(fig_batch)==3){
        assemble_plot <- wrap_plots(fig_batch,ncol=1)
        pdf_name <- paste0("Boxplot_Batch_",length(list.files(pattern="Boxplot_Batch*"))+1,".pdf")
        pdf(pdf_name,width=9.5,height=11)
        print(assemble_plot)
        dev.off()
        fig_batch <- list()
      }
    }
  }
}

if(length(fig_batch)>0){
  assemble_plot <- wrap_plots(fig_batch,ncol=1)
  pdf("Boxplot_Batch_remain.pdf",width=9.5,height=11)
  print(assemble_plot)
  dev.off()
}

stat_res_final <- stat_res_all %>% mutate(kruskal_p_BH=p.adjust(kruskal_p,method="BH"),cell_level_p_BH=p.adjust(cell_level_p,method="BH"),gene_group_level_p_BH=p.adjust(gene_group_level_p,method="BH"))
write_csv(stat_res_final,"./scoreBox/Pathway_CellType_StatResult.csv")

plot_list$Necroptosis___NK
plot_list$`TNFA_signaling_via_NFKB___Non-naive_CD8_T`
plot_list$`TNFA_signaling_via_NFKB___NK`
plot_list$`TNFA_signaling_via_NFKB___Plasma`
plot_list$Interferon_stimulated_genes___Plasma



p1 <- plot_list[["Necroptosis___NK"]]
p2 <- plot_list[["TNFA_signaling_via_NFKB___Non-naive_CD8_T"]]
p3 <- plot_list[["TNFA_signaling_via_NFKB___NK"]]
p4 <- plot_list[["TNFA_signaling_via_NFKB___Plasma"]]
p5 <- plot_list[["Interferon_stimulated_genes___Plasma"]]

style_upper <- function(p){
  p + 
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    )
}

p1 <- style_upper(p1)
p2 <- style_upper(p2)
p3 <- style_upper(p3)
p4 <- style_upper(p4)

composite <- (p1 / p2 / p3 / p4 / p5) +
  plot_layout(ncol = 1, heights = rep(1,5)) +
  plot_annotation(tag_levels = NULL) &
  theme(plot.margin = margin(0,0,0,0))

ggsave("./scoreBox/selected_combined_panel.pdf", composite, width = 8.5, height = 13.5)


