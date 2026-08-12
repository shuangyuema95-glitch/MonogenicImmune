library(ggplot2)
library(patchwork)
library(viridis)
library(tidyverse)


cond_vec <- c("Arachidonic acid metabolism","Ca2+_flux-PLC","Cytoskeleton_and_small_GTPase",
              "Endolysosomal_nucleic_acid_sensing","HC","Immune_metabolic",
              "Inborn_errors_of_cell_death","Inflammasome_IL-1β","JIA","KD",
              "Negative_regulation_of_IFN-I","NFKB pathway","Osteoclast function",
              "Protein_homeostasis","SLE","TBK1_IRF3","Uncategorized")

cond_order <- c("HC","KD","JIA","SLE","Uncategorized","Protein_homeostasis",
                "TBK1_IRF3","Cytoskeleton_and_small_GTPase","Negative_regulation_of_IFN-I","Inborn_errors_of_cell_death",
                "Immune_metabolic","Inflammasome_IL-1β","Arachidonic acid metabolism",
                "Endolysosomal_nucleic_acid_sensing","Ca2+_flux-PLC","Osteoclast function","NFKB pathway")

cond_rename_tbl <- tibble(
  full = cond_order,
  short = c("HC","KD","JIA","SLE","UNC","PHS",
            "TIR","CSG","NRI","IECD","IMM","IIL","AAM",
            "ENS","CFP","OST","NFB"))

#' CytoTRACE2 Violin Plot Wrapper
#' Two‑mode violin plot for CytoTRACE2 score.
#' mode = "condimodel": loop over conditions, x‑axis = cell subclusters, fill by cell_color_map.
#' mode = "cellmodel": loop over cell subclusters, x‑axis = conditions, unified fill by cond_violin_fill.
#'
#' @param Mono_sub Seurat object list, contains subMeta slot with barcode and Level2 metadata
#' @param Cyto_mono data.frame, CytoTRACE2 output metadata, first column is barcode
#' @param mode character, one of c("cellmodel","condimodel")
#' @param cond_vec character vector, condition names to iterate for plotting
#' @param cell_order character vector, factor levels for Level2 cell subclusters
#' @param cell_color_map named character vector, color mapping for Level2 subclusters
#' @param cond_order character vector, factor levels for condition
#' @param cond_violin_fill character, fill color for violin under cellmodel mode
#'
#' @return list, list(fig_list = list of individual ggplot objects, full_plot = wrapped patchwork object)
#' @export


cyto_trace_plot <- function(sub,
                            cyto,
                            mode = c("condimodel_violin",
                                     "cellmodel_violin",
                                     "condimodel_boxonly",
                                     "cellmodel_boxonly"),
                            cond_vec,
                            cell_order,
                            cell_color_map,
                            cond_order,
                            cond_violin_fill = "#0B71AB",
                            cond_box_color = "#0B71AB",
                            p_adj_method = "BH"){
  mode <- match.arg(mode)
  meta_df <- as.data.frame(sub$subMeta)
  meta_df$barcode <- rownames(meta_df)
  meta_df$Level2 <- factor(meta_df$Level2, levels = cell_order)
  cyto_df <- cyto
  colnames(cyto_df)[1] <- "barcode"
  plt_all <- merge(meta_df[,c("barcode","Level2"),drop=FALSE],cyto_df[,c("barcode","condition","CytoTRACE2_Score"),drop=FALSE],by="barcode")
  plt_all$condition <- factor(plt_all$condition,levels=cond_order)
  
  cond_color_vec <- dplyr::tibble(condition = cond_order) %>%
    dplyr::mutate(color = dplyr::case_when(
      condition == "HC" ~ "#4375A0",
      condition %in% c("KD","JIA","SLE") ~ "#B18749",
      TRUE ~ "#912623"
    )) %>%
    dplyr::pull(color, name = condition)
  
  cond_kw_raw <- list()
  for(cond in cond_vec){
    dtmp <- plt_all[plt_all$condition == cond, , drop = FALSE]
    if(nrow(dtmp) < 3 || length(unique(dtmp$Level2)) < 2){
      cond_kw_raw[[cond]] <- tibble(group_id = cond, p_value = NA_real_)
      next
    }
    kt <- kruskal.test(CytoTRACE2_Score ~ Level2, data = dtmp)
    cond_kw_raw[[cond]] <- tibble(group_id = cond, p_value = kt$p.value)
  }
  cond_kw_df <- bind_rows(cond_kw_raw) %>% mutate(p_adj = p.adjust(p_value, method = p_adj_method),test_type = "within_condition_across_cells") %>% rename(condition = group_id)
  
  cell_kw_raw <- list()
  for(ct in cell_order){
    dtmp <- plt_all[plt_all$Level2 == ct, , drop = FALSE]
    if(nrow(dtmp) < 3 || length(unique(dtmp$condition)) < 2){
      cell_kw_raw[[ct]] <- tibble(group_id = ct, p_value = NA_real_)
      next
    }
    kt <- kruskal.test(CytoTRACE2_Score ~ condition, data = dtmp)
    cell_kw_raw[[ct]] <- tibble(group_id = ct, p_value = kt$p.value)
  }
  cell_kw_df <- bind_rows(cell_kw_raw) %>% mutate(p_adj = p.adjust(p_value, method = p_adj_method),test_type = "within_cell_across_conditions") %>% rename(cell_type = group_id)
  
  if(mode %in% c("condimodel_violin","condimodel_boxonly")){
    cell_label_map <- purrr::set_names(stringr::str_extract(cell_order, "^[^_]+"),nm = cell_order)
    out_list <- list()
    for(cond in cond_vec){
      cyto_sub <- cyto_df[cyto_df$condition == cond, , drop = FALSE]
      cell_vec <- cyto_sub$barcode
      if(length(cell_vec)==0){message(sprintf("Skip: %s, zero cells",cond));next}
      plt_df <- merge(meta_df[meta_df$barcode %in% cell_vec,c("barcode","Level2"),drop=FALSE],cyto_sub[,c("barcode","CytoTRACE2_Score"),drop=FALSE],by="barcode")
      if(nrow(plt_df)<1){message(sprintf("Skip: %s, no overlapping cells",cond));next}
      if(mode == "condimodel_violin"){
        p <- ggplot(plt_df,aes(x=Level2,y=CytoTRACE2_Score,fill=Level2))+
          geom_violin(trim=FALSE,scale="width",width=0.32,size=0,color=NA)+
          geom_boxplot(width=0.35,color="black",outlier.shape=NA,fill=NA,fatten=0.3)+
          labs(title=cond,y="CytoTRACE2_Score",x="")+
          scale_fill_manual(values=cell_color_map)+
          scale_x_discrete(labels = cell_label_map)+
          theme_classic()+
          theme(plot.title=element_text(hjust=0.5,colour="black"),legend.position="none",axis.ticks=element_line(color="black"),axis.text.x=element_text(colour="black",angle=45,hjust=1),axis.text.y=element_text(colour="black"),axis.title=element_text(colour="black"))
      } else {
        p <- ggplot(plt_df,aes(x=Level2,y=CytoTRACE2_Score,color=Level2))+
          geom_boxplot(fill=NA,width=0.55,linewidth=0.6,outlier.shape=NA)+
          scale_color_manual(values=cell_color_map)+
          scale_x_discrete(labels = cell_label_map)+
          labs(title=cond,y="CytoTRACE2_Score",x="")+
          theme_classic()+
          theme(plot.title=element_text(hjust=0.5,colour="black"),legend.position="none",axis.ticks=element_line(color="black"),axis.text.x=element_text(colour="black",angle=45,hjust=1),axis.text.y=element_text(colour="black"),axis.title=element_text(colour="black"))
      }
      out_list[[cond]] <- p
    }
    full_p <- patchwork::wrap_plots(out_list)
    return(list(fig_list = out_list,full_plot = full_p,cond_kw_df = cond_kw_df,cell_kw_df = cell_kw_df,cell_cond_mean_list=list()))
  }
  
  if(mode %in% c("cellmodel_violin","cellmodel_boxonly")){
    out_list <- list()
    cell_cond_mean_list <- list()
    for(ct in cell_order){
      sub_df <- plt_all[plt_all$Level2 == ct,,drop=FALSE]
      keep_cond <- names(table(sub_df$condition))[table(sub_df$condition)>0]
      sub_df <- sub_df[sub_df$condition %in% keep_cond,,drop=FALSE]
      if(nrow(sub_df)<1){message(sprintf("Skip %s: zero cells",ct));next}
      sub_df$condition <- factor(sub_df$condition,levels=cond_order[cond_order %in% keep_cond])
      
      hc_median_percell <- sub_df %>% dplyr::filter(condition=="HC") %>% dplyr::pull(CytoTRACE2_Score) %>% median(na.rm=TRUE)
      
      percell_df <- sub_df %>%
        dplyr::group_by(condition) %>%
        dplyr::summarise(group_median = median(CytoTRACE2_Score,na.rm=TRUE),.groups="drop") %>%
        dplyr::filter(condition=="HC" | group_median < hc_median_percell)
      cell_cond_mean_list[[ct]] <- percell_df
      
      current_lv <- levels(sub_df$condition)
      label_short <- cond_rename_tbl$short[match(current_lv, cond_rename_tbl$full)]
      
      if(mode == "cellmodel_violin"){
        p <- ggplot(sub_df,aes(x=condition,y=CytoTRACE2_Score,fill=condition))+
          geom_hline(yintercept = hc_median_percell, color = "gray60", linetype = 2, linewidth = 0.4)+
          geom_violin(trim=FALSE,scale="width",width=0.32,size=0,color=NA)+
          geom_boxplot(width=0.35,color="black",outlier.shape=NA,fill=NA,fatten=0.3)+
          scale_fill_manual(values = cond_color_vec)+
          scale_x_discrete(labels = label_short)+
          labs(title=ct,y="CytoTRACE2_Score",x="")+
          theme_classic()+
          theme(plot.title=element_text(hjust=0.5,colour="black"),legend.position="none",axis.ticks=element_line(color="black"),axis.text.x=element_text(colour="black",angle=45,hjust=1),axis.text.y=element_text(colour="black"),axis.title=element_text(colour="black"))
      } else {
        p <- ggplot(sub_df,aes(x=condition,y=CytoTRACE2_Score,color=condition))+
          geom_hline(yintercept = hc_median_percell, color = "gray60", linetype = 2, linewidth = 0.4)+
          geom_boxplot(fill=NA,width=0.55,linewidth=0.6,outlier.shape=NA)+
          scale_color_manual(values = cond_color_vec)+
          scale_x_discrete(labels = label_short)+
          labs(title=ct,y="CytoTRACE2_Score",x="")+
          theme_classic()+
          theme(plot.title=element_text(hjust=0.5,colour="black"),legend.position="none",axis.ticks=element_line(color="black"),axis.text.x=element_text(colour="black",angle=45,hjust=1),axis.text.y=element_text(colour="black"),axis.title=element_text(colour="black"))
      }
      out_list[[ct]] <- p
    }
    full_p <- patchwork::wrap_plots(out_list)
    return(list(fig_list = out_list,full_plot = full_p,cond_kw_df = cond_kw_df,cell_kw_df = cell_kw_df,cell_cond_mean_list=cell_cond_mean_list))
  }
}


##(1) Monocyte 
Cyto_mono<-read.csv("E:\\AID cohort\\code\\Cellular_trajectory\\monocyte_cytotrace2_metadata.csv")
dim(Cyto_mono)
Mono_sub<-readRDS("F:\\AID_results_backup\\subseurat_rds\\Mono_subresult.RDS")
names(Mono_sub)
cell_order <- c("Mono01_Classical_FOSB","Mono02_Inflammatory","Mono03_Nonclassical",
                "Mono04_Classical_RETN","Mono05_Classical_S100A9","Mono06_Classical_TLN1","Mono07_IFNresponse")

cell_color_map <- c("Mono01_Classical_FOSB"="#847AB3","Mono02_Inflammatory"="#003D81","Mono03_Nonclassical"="#C6D199",
                    "Mono04_Classical_RETN"="#FFEFC1","Mono05_Classical_S100A9"="#F6E1EE","Mono06_Classical_TLN1"="#D1929B","Mono07_IFNresponse"="#62a7bd")

res_cell_box <- cyto_trace_plot(sub = Mono_sub,
                                cyto = Cyto_mono,
                                mode = "cellmodel_boxonly",
                                cond_vec = cond_vec,
                                cell_order = cell_order,
                                cell_color_map = cell_color_map,
                                cond_order = cond_order,
                                cond_box_color = "#0B71AB")

res_cell_box$cond_kw_df
res_cell_box$cell_kw_df%>%filter( p_adj<0.01)# 1 2 3 4 5 6 (0.05,same)
res_cell_box$fig_list$Mono01_Classical_FOSB;res_cell_box$cell_cond_mean_list$Mono01_Classical_FOSB
res_cell_box$fig_list$Mono02_Inflammatory;res_cell_box$cell_cond_mean_list$Mono02_Inflammatory
res_cell_box$fig_list$Mono03_Nonclassical;res_cell_box$cell_cond_mean_list$Mono03_Nonclassical
res_cell_box$fig_list$Mono04_Classical_RETN;res_cell_box$cell_cond_mean_list$Mono04_Classical_RETN
res_cell_box$fig_list$Mono05_Classical_S100A9;res_cell_box$cell_cond_mean_list$Mono05_Classical_S100A9
res_cell_box$fig_list$Mono06_Classical_TLN1;res_cell_box$cell_cond_mean_list$Mono06_Classical_TLN1

P1_cyto<-res_cell_box$fig_list$Mono01_Classical_FOSB+theme(
                                                           axis.title.y = element_blank(),
                                                           plot.title = element_blank())
P3_cyto<-res_cell_box$fig_list$Mono03_Nonclassical+theme(axis.title.y = element_blank(),
                                                         plot.title = element_blank())
Cyto<-P1_cyto+P3_cyto
Cyto


ggsave(Cyto,file="E:\\AID cohort\\code\\Cellular_trajectory\\Mono0103_Cyto.pdf",
       width=8.09,height = 2.48)

##################
###scatter for specific
library(tidyverse)
library(patchwork)

cond_order <- c("HC","KD","JIA","SLE","Uncategorized","Protein_homeostasis",
                "TBK1_IRF3","Cytoskeleton_and_small_GTPase","Negative_regulation_of_IFN-I","Inborn_errors_of_cell_death",
                "Immune_metabolic","Inflammasome_IL-1β","Arachidonic acid metabolism",
                "Endolysosomal_nucleic_acid_sensing","Ca2+_flux-PLC","Osteoclast function","NFKB pathway")
cond_rename_tbl <- tibble(full=cond_order,short=c("HC","KD","JIA","SLE","UNC","PHS","TIR","CSG","NRI","IECD","IMM","IIL","AAM","ENS","CFP","OST","NFB"))

#(1) different color
gene_loop_patch <- function(seu_obj,target_ct,gene_vec,gene_color_vec,cond_order,cond_rename_df,stat_type=c("mean","median"),plot_mode=c("bar","line"),add_smooth=FALSE,bar_width=0.7,point_size=2,line_width=0.6){
  stat_type<-match.arg(stat_type);plot_mode<-match.arg(plot_mode)
  seu_work<-seu_obj;seu_work@meta.data$condition<-factor(seu_work@meta.data$condition,levels=cond_order)
  meta_df<-seu_work@meta.data%>%rownames_to_column("barcode")%>%filter(Level2==target_ct)
  expr_mat<-FetchData(seu_work,vars=gene_vec)%>%rownames_to_column("barcode")
  comb_df<-inner_join(meta_df,expr_mat,by="barcode")
  sum_fun<-switch(stat_type,mean=~mean(.x,na.rm=TRUE),median=~median(.x,na.rm=TRUE))
  plot_df<-comb_df%>%group_by(condition)%>%summarise(across(all_of(gene_vec),sum_fun),.groups="drop")%>%pivot_longer(cols=all_of(gene_vec),names_to="gene",values_to="expr_val")%>%mutate(condition=factor(condition,levels=cond_order))%>%left_join(cond_rename_df,by=c("condition"="full"))
  plot_df$short<-factor(plot_df$short,levels=cond_rename_df$short);plot_df$x_pos<-as.integer(plot_df$short)
  p_list<-list()
  for(g in gene_vec){
    sub_df<-plot_df%>%filter(gene==g)
    hc_ref<-sub_df%>%filter(condition=="HC")%>%pull(expr_val)
    p<-ggplot(sub_df,aes(x=x_pos,y=expr_val))+
      geom_hline(yintercept=hc_ref,linetype="dashed",color="#666666",linewidth=0.4)+
      labs(x=NULL,y=g)+
      scale_x_continuous(breaks=seq_along(cond_rename_df$short),labels=cond_rename_df$short)+
      theme_classic()+
      theme(panel.grid=element_blank(),axis.text=element_text(colour="black"),axis.title=element_text(colour="black"),axis.text.x=element_text(angle=45,hjust=1))
    if(plot_mode=="bar"){
      p<-p+geom_col(fill=gene_color_vec[[g]],width=bar_width)
    }else{
      p<-p+geom_line(aes(group=1),color=gene_color_vec[[g]],linewidth=line_width)+geom_point(color=gene_color_vec[[g]],size=point_size)
      if(isTRUE(add_smooth)){p<-p+geom_smooth(aes(group=1),se=FALSE,method="loess",span=0.8,color=gene_color_vec[[g]],linewidth=line_width)}
    }
    p_list[[g]]<-p
  }
  n_p<-length(p_list)
  for(i in seq_along(p_list)){if(i<n_p){p_list[[i]]<-p_list[[i]]+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())}}
  p_all<-wrap_plots(p_list,ncol=1,heights=rep(1,n_p))
  return(p_all)
}

gene_color_vec<-c(CD86="#E41A1C",IL1B="#377EB8",IRF1="#4DAF4A",KYNU="#984EA3")
p_bar_all<-gene_loop_patch(seu_obj=Mono_seu,target_ct="Mono02_Inflammatory",gene_vec=c("CD86","IL1B","IRF1","KYNU"),gene_color_vec=gene_color_vec,cond_order=cond_order,cond_rename_df=cond_rename_tbl,stat_type="mean",plot_mode="bar",bar_width=0.7)
p_line_all<-gene_loop_patch(seu_obj=Mono_seu,target_ct="Mono02_Inflammatory",gene_vec=c("CD86","IL1B","IRF1","KYNU"),gene_color_vec=gene_color_vec,cond_order=cond_order,cond_rename_df=cond_rename_tbl,stat_type="mean",plot_mode="line",add_smooth=FALSE,point_size=2,line_width=0.6)

p_bar_all
p_line_all

#same color
library(tidyverse)
library(patchwork)

#' gene_loop_patch
#' @description Generate vertically‑aligned multi‑panel bar or line plots for multiple genes.
#' Compute mean/median expression per condition, add dashed reference line for HC group.
#' @param seu_obj Seurat object
#' @param target_ct Character, Level2 cell type for subsetting cells
#' @param gene_vec Character vector, target gene symbols to plot

gene_loop_patch <- function(seu_obj,target_ct,gene_vec,base_color="#912623",
                            Lower_color="#f8c8c8",cond_order,cond_rename_df,
                            stat_type=c("mean","median"),plot_mode=c("bar","line"),
                            add_smooth=FALSE,remove_original_line_when_smooth=FALSE,
                            bar_width=0.7,point_size=0.89,line_width=0.4,axis_text_size=6){
  stat_type<-match.arg(stat_type);plot_mode<-match.arg(plot_mode)
  seu_work<-seu_obj;seu_work@meta.data$condition<-factor(seu_work@meta.data$condition,levels=cond_order)
  meta_df<-seu_work@meta.data%>%rownames_to_column("barcode")%>%filter(Level2==target_ct)
  expr_mat<-FetchData(seu_work,vars=gene_vec)%>%rownames_to_column("barcode")
  comb_df<-inner_join(meta_df,expr_mat,by="barcode")
  sum_fun<-switch(stat_type,mean=~mean(.x,na.rm=TRUE),median=~median(.x,na.rm=TRUE))
  plot_df<-comb_df%>%group_by(condition)%>%summarise(across(all_of(gene_vec),sum_fun),.groups="drop")%>%pivot_longer(cols=all_of(gene_vec),names_to="gene",values_to="expr_val")%>%mutate(condition=factor(condition,levels=cond_order))%>%left_join(cond_rename_df,by=c("condition"="full"))
  plot_df$short<-factor(plot_df$short,levels=cond_rename_df$short);plot_df$x_pos<-as.integer(plot_df$short)
  p_list<-list()
  for(g in gene_vec){
    sub_df<-plot_df%>%filter(gene==g)
    hc_ref<-sub_df%>%filter(condition=="HC")%>%pull(expr_val)
    p<-ggplot(sub_df,aes(x=x_pos,y=expr_val))+
      geom_hline(yintercept=hc_ref,linetype="dashed",color="#666666",linewidth=0.4)+
      labs(x=NULL,y=g)+
      scale_x_continuous(breaks=seq_along(cond_rename_df$short),labels=cond_rename_df$short)+
      theme_classic()+
      theme(panel.grid=element_blank(),
            axis.text=element_text(colour="black",size=axis_text_size),
            axis.title=element_text(colour="black"),
            axis.text.x=element_text(angle=45,hjust=1),
            legend.position="none")
    if(plot_mode=="bar"){
      p<-p+geom_col(aes(fill=expr_val),width=bar_width)+
        scale_fill_gradient(low=Lower_color,high=base_color)
    }else{
      if(!(add_smooth && remove_original_line_when_smooth)){
        p<-p+geom_line(aes(group=1),color=base_color,linewidth=line_width)
      }
      p<-p+geom_point(color=base_color,size=point_size)
      if(isTRUE(add_smooth)){p<-p+geom_smooth(aes(group=1),se=FALSE,method="loess",span=0.8,color=base_color,linewidth=line_width)}
    }
    p_list[[g]]<-p
  }
  n_p<-length(p_list)
  #for(i in seq_along(p_list)){if(i<n_p){p_list[[i]]<-p_list[[i]]+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())}}
  p_all<-wrap_plots(p_list,nrow=1,heights=rep(1,n_p))
  return(p_all)
}

# call
p_line_all1<-gene_loop_patch(seu_obj=Mono_seu,target_ct="Mono01_Classical_FOSB",gene_vec=c("CD86","KYNU","IL1B"),
                             base_color="#8C2522",cond_order=cond_order,
                             cond_rename_df=cond_rename_tbl,stat_type="mean",
                             plot_mode="line",add_smooth=TRUE,remove_original_line_when_smooth=TRUE,point_size=0.1,line_width=0.25,axis_text_size=6)

p_line_all1

p_line_all2<-gene_loop_patch(seu_obj=Mono_seu,target_ct="Mono03_Nonclassical",gene_vec=c("CD86","KYNU","IL15"),
                             base_color="#8C2522",cond_order=cond_order,
                             cond_rename_df=cond_rename_tbl,stat_type="mean",
                             plot_mode="line",add_smooth=TRUE,remove_original_line_when_smooth=TRUE,point_size=0.1,line_width=0.45,axis_text_size=6)
p_line_all2

expMacroPro<-p_line_all1/p_line_all2

expMacroPro
ggsave(expMacroPro,file="E:\\AID cohort\\code\\Cellular_trajectory\\expMacroPro.pdf",
       width=8.09,height = 2.48)

#选择IL15 支撑文献:IL-15 functions as a potent autocrine regulator of macrophage proinflammatory cytokine production: evidence for differential receptor subunit utilization associated with stimulation or inhibition



# p_bar_all2<-gene_loop_patch(seu_obj=Mono_seu,
#                             target_ct="Mono01_Classical_FOSB",gene_vec=c("CD86","IL1B","KYNU"),
#                             base_color="#E41A1C",Lower_color = "#FFEFC1",cond_order=cond_order,
#                             cond_rename_df=cond_rename_tbl,stat_type="mean",plot_mode="bar",bar_width=0.65,axis_text_size=6)
# 
# p_bar_all2






##(2) DC
# Cyto_DC<-read.csv("E:\\AID cohort\\code\\Cellular_trajectory\\DC_cytotrace2_metadata.csv")
# dim(Cyto_DC)
# DC_sub<-readRDS("F:\\AID_results_backup\\subseurat_rds\\DC_subresult.RDS")
# names(DC_sub)
# 
# cell_order<-sort(unique(DC_sub$subMeta$Level2))
# colors=c("#847AB3","#003D81","#C6D199","#FFEFC1","#F6E1EE","#D1929B","#62a7bd","#0B71AB","#C2D4F2")
# 
# cell_color_map<-colors
# names(cell_color_map)<-cell_order
# cell_color_map
# 
# 
# res_cell_box2 <- cyto_trace_plot(sub = DC_sub,
#                                 cyto = Cyto_DC,
#                                 mode = "cellmodel_boxonly",
#                                 cond_vec = cond_vec,
#                                 cell_order = cell_order,
#                                 cell_color_map = cell_color_map,
#                                 cond_order = cond_order,
#                                 cond_box_color = "#0B71AB")


# res_condi_box <- cyto_trace_plot(sub = DC_sub,
#                                  cyto = Cyto_DC,
#                                  mode = "condimodel_boxonly",
#                                  cond_vec = cond_vec,
#                                  cell_order = cell_order,
#                                  cell_color_map = cell_color_map)


# res_cell_box2$cond_kw_df
# res_cell_box2$cell_kw_df%>%filter( p_adj<0.01)# 1 2 3 4 5 6 8
# res_cell_box2$fig_list$DC01_cDC2_CD1C;res_cell_box2$cell_cond_mean_list$DC01_cDC2_CD1C
# res_cell_box2$fig_list$DC02_DC_TRIM33;res_cell_box2$cell_cond_mean_list$DC02_DC_TRIM33
# res_cell_box2$fig_list$DC03_aDC_CCR7;res_cell_box2$cell_cond_mean_list$DC03_aDC_CCR7
# res_cell_box2$fig_list$DC04_DC_CCL5;res_cell_box2$cell_cond_mean_list$DC04_DC_CCL5
# res_cell_box2$fig_list$DC05_cDC2_LY96;res_cell_box2$cell_cond_mean_list$DC05_cDC2_LY96
# res_cell_box2$fig_list$DC06_cDC1_XCR1;res_cell_box2$cell_cond_mean_list$DC06_cDC1_XCR1#use this
# res_cell_box2$fig_list$DC08_DC_C1QTNF4;res_cell_box2$cell_cond_mean_list$DC08_DC_C1QTNF4


###focused on lowset CV value in Figure 2A (broadly involved in DEGs identification)
##(3) pDC
##(4) LDG
##(5) Plasma
##(6) UTC


##################
###violin to check anti-inflammation or pro-inflammation gene expressions
library(Seurat)
library(data.table)
library(tidyverse)
Mono_seu<-readRDS("F:\\AID_results_backup\\subseurat_rds\\Mono_level2_seu.rds")
Idents(Mono_seu)<-"level2_cluster"
Mono_seu<-RenameIdents(Mono_seu,
                       "0"="Mono01_Classical_FOSB",
                       "1"="Mono02_Inflammatory",
                       "2"="Mono03_Nonclassical",
                       "3"="Mono04_Classical_RETN",
                       "4"="Mono05_Classical_S100A9",
                       "5"="Mono06_Classical_TLN1",
                       "6"="Mono01_Classical_FOSB",
                       "7"="Mono07_IFNresponse",
                       "8"="Mono07_IFNresponse",
                       "9"="Mono02_Inflammatory")
Mono_seu@meta.data$Level2 <- Idents(Mono_seu)

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
meta<-Mono_seu@meta.data
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
table(meta$condition)
Mono_seu@meta.data<-meta;rm(meta);gc()

ImmuneList<- as.data.frame(fread("E:\\AID cohort\\data\\0725_FinalList.txt"))
# ImmuneList<-ImmuneList[ImmuneList$name%in%c("Cytokine_receptors_antiinflammatory",
#                                             "Cytokine_receptors_proinflammatory","IL10_antiinflammatory_signaling_pathway",
#                                             "M1_macrophage","M2_macrophage"),c("name","gene")]

ImmuneList<-ImmuneList[ImmuneList$name%in%c("M1_macrophage","M2_macrophage"),c("name","gene")]


genelist<-list()
for(i in 1:nrow(ImmuneList)){
  genes<-trimws(unlist(strsplit(ImmuneList[i,2],split = ",")))
  genes<-rownames(Mono_seu)[na.omit(match(genes,rownames(Mono_seu)))]
  genelist[[i]]<-genes
}
names(genelist)<-ImmuneList$name
#names(genelist)<-c("anti_cytok","pro_cytok","IL10_anti","M1_macro","M2_macro")
names(genelist)<-c("M1_macro","M2_macro")

sigcyto_cluster<-cell_order <- c("Mono01_Classical_FOSB","Mono02_Inflammatory","Mono03_Nonclassical",
                                 "Mono04_Classical_RETN","Mono05_Classical_S100A9","Mono06_Classical_TLN1")


cond_order <- c("HC","KD","JIA","SLE","Uncategorized","Protein_homeostasis",
                "TBK1_IRF3","Cytoskeleton_and_small_GTPase","Negative_regulation_of_IFN-I","Inborn_errors_of_cell_death",
                "Immune_metabolic","Inflammasome_IL-1β","Arachidonic acid metabolism",
                "Endolysosomal_nucleic_acid_sensing","Ca2+_flux-PLC","Osteoclast function","NFKB pathway")

cond_rename_tbl <- tibble(
  full = cond_order,
  short = c("HC","KD","JIA","SLE","UNC","PHS",
            "TIR","CSG","NRI","IECD","IMM","IIL","AAM",
            "ENS","CFP","OST","NFB"))

cond_color <- dplyr::case_when(
  cond_order == "HC" ~ "#4375A0",
  cond_order %in% c("KD","JIA","SLE") ~ "#B18749",
  TRUE ~ "#912623"
)
names(cond_color) <- cond_order
table(Mono_seu$condition)


viogene_mod <- function(seu_obj,
                        intercell,
                        gene_df,
                        cond_color_vec,
                        cond_rename){
  
  data_exp <- GetAssayData(seu_obj, assay = "RNA", slot = "data")
  
  meta_raw <- seu_obj@meta.data %>%
    rownames_to_column("cell_id") %>%
    filter(Level2 == intercell, !is.na(condition))
  
  meta_raw <- meta_raw %>%
    mutate(cond_fac = factor(condition, levels = cond_rename$full))
  
  plot_list <- list()
  
  for (k in seq_len(nrow(gene_df))) {
    gene <- gene_df$gene[k]
    setname <- gene_df$set_name[k]
    
    gene_idx <- match(gene, rownames(data_exp))
    if(is.na(gene_idx)){
      message(sprintf("[%s] gene %s not found, skip", intercell, gene))
      next
    }
    
    cell_vec <- meta_raw$cell_id
    expr_vec <- as.numeric(data_exp[gene_idx, cell_vec])
    
    if(all(expr_vec == 0, na.rm = TRUE)){
      message(sprintf("[%s] gene %s all zero, skip", intercell, gene))
      next
    }
    
    plot_df <- meta_raw %>%
      mutate(value = expr_vec)
    
    p <- ggplot(plot_df, aes(x = cond_fac, y = value)) +
      geom_violin(aes(fill = cond_fac), scale = "width", adjust = 1, trim = TRUE) +
      scale_fill_manual(values = cond_color_vec) +
      scale_x_discrete(limits = cond_rename$full, labels = cond_rename$short) +
      labs(title = gene, y = setname, x = NULL) +
      theme_bw() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.title = element_text(face = "italic", size = 12, hjust = 0.5),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA)
      )
    plot_list[[length(plot_list)+1]] <- p
  }
  return(plot_list)
}


gene_tbl <- imap_dfr(genelist, ~tibble(gene = .x, set_name = .y))
gene_tbl <- gene_tbl %>% filter(gene %in% rownames(Mono_seu@assays$RNA@data))

mono_violin_res <- list()
for(one_cell in sigcyto_cluster){
  message("\n==== Start process Level2: ", one_cell, " ====")
  
  if(!one_cell %in% unique(Mono_seu$Level2)){
    message("Warning: ", one_cell, " not present in object, skip")
    next
  }
  
  gg_plots <- viogene_mod(
    seu_obj = Mono_seu,
    intercell = one_cell,
    gene_df = gene_tbl,
    cond_color_vec = cond_color,
    cond_rename = cond_rename_tbl
  )
  
  mono_violin_res[[one_cell]] <- list(gg_plots = gg_plots)
  message("---- ", one_cell, " valid plot count: ", length(gg_plots))
}

##saving
library(patchwork)
n_col <- 3
n_row <- 5
per_page <- n_col * n_row

for(ct_name in names(mono_violin_res)){
  plist <- mono_violin_res[[ct_name]]$gg_plots
  if(length(plist) == 0){
    message(sprintf("\n%s : no valid plots, skip pdf output", ct_name))
    next
  }
  pdf_name <- paste0(gsub("[| :]", "_", ct_name), "_violin.pdf")
  message(sprintf("\nWrite pdf: %s, total plots: %d", pdf_name, length(plist)))
  pdf(pdf_name, width = 16, height = 12)
  split_idx <- split(seq_along(plist), ceiling(seq_along(plist)/per_page))
  for(page_i in names(split_idx)){
    sub_plist <- plist[split_idx[[page_i]]]
    page_comb <- wrap_plots(sub_plist, ncol = n_col, nrow = n_row)
    print(page_comb)
  }
  dev.off()
}


##################
###Bubble plot to check anti-inflammation or pro-inflammation gene expressions
library(tidyverse)
library(RColorBrewer)


cond_order <- c("HC","KD","JIA","SLE","Uncategorized","Protein_homeostasis",
                "TBK1_IRF3","Cytoskeleton_and_small_GTPase","Negative_regulation_of_IFN-I","Inborn_errors_of_cell_death",
                "Immune_metabolic","Inflammasome_IL-1β","Arachidonic acid metabolism",
                "Endolysosomal_nucleic_acid_sensing","Ca2+_flux-PLC","Osteoclast function","NFKB pathway")

cond_rename_tbl <- tibble(
  full = cond_order,
  short = c("HC","KD","JIA","SLE","UNC","PHS",
            "TIR","CSG","NRI","IECD","IMM","IIL","AAM",
            "ENS","CFP","OST","NFB"))



gene_mean_square_bubble <- function(seu_obj,
                                    target_ct,
                                    gene_vec,
                                    cond_order,
                                    cond_rename_df,
                                    pal_gradient,
                                    stat_type = c("mean","median")){
  
  stat_type <- match.arg(stat_type)
  seu_work <- seu_obj
  seu_work@meta.data$condition <- factor(seu_work@meta.data$condition, levels = cond_order)
  
  meta_df <- seu_work@meta.data %>%
    rownames_to_column("barcode") %>%
    filter(Level2 == target_ct)
  
  expr_mat <- FetchData(seu_work, vars = gene_vec) %>%
    rownames_to_column("barcode")
  
  comb_df <- inner_join(meta_df, expr_mat, by = "barcode")
  
  sum_fun <- switch(stat_type,
                    mean = ~mean(.x, na.rm = TRUE),
                    median = ~median(.x, na.rm = TRUE))
  
  plot_df <- comb_df %>%
    group_by(condition) %>%
    summarise(across(all_of(gene_vec), sum_fun), .groups = "drop") %>%
    pivot_longer(cols = all_of(gene_vec), names_to = "gene", values_to = "expr_val") %>%
    mutate(
      condition = factor(condition, levels = cond_order),
      gene = factor(gene, levels = gene_vec)
    ) %>%
    left_join(cond_rename_df, by = c("condition" = "full"))
  
  plot_df$short <- factor(plot_df$short, levels = cond_rename_df$short)
  
  p <- ggplot(plot_df, aes(x = short, y = gene)) +
    geom_point(aes(size = expr_val, fill = expr_val), shape = 22, stroke = 0) +
    scale_fill_gradientn(colors = pal_gradient) +
    scale_size_continuous(range = c(2,10)) +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(
      panel.grid = element_line(color = "grey80"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      axis.line = element_blank(),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
  return(p)
}

col_raw <- c("#FFFFFF","#F7F7F8","#EFEDEE","#F1DFD8","#E6A77C","#D17C42","#B3563C","#97332D","#832528")
col_100 <- colorRampPalette(col_raw)(100)

col_raw <- c(
  "#BFBFBF","#D6D6C8","#E9F6B1","#B7DDAF",
  "#7AC0BB","#45B3C2","#2195C0","#1E7FB7",
  "#215FA8","#154092","#0C2C84")
col_100 <- colorRampPalette(col_raw)(100)
col_100

# mean
p_mean <- gene_mean_square_bubble(
  seu_obj = Mono_seu,
  target_ct = "Mono01_Classical_FOSB",
  gene_vec = c("CD86","IL1B","IRF1","KYNU"),
  cond_order = cond_order,
  cond_rename_df = cond_rename_tbl,
  pal_gradient = col_100,
  stat_type = "mean"
)

