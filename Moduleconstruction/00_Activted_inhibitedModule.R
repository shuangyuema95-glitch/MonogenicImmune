#======================================================================
# Disease Module NMF — Separate Activation/Inhibition Networks
# Input:  sigLeme_widebeta.RDS (174 factors x 16 diseases, LMEM betas)
# Output: 0902_K12NMF.Rdata, 0902_NMF.Rdata (full workspace)
# Steps:  1) Filter+celltype  2) NMF(K=12)  3) Spearman+specificity
#         4) Sensitivity top_n (act/inh separate, BH internal)
#         5) Build edge lists (peak top_n, act/inh separate)
#         6) Module activity + bidding (act/inh)
#         7) Resolve conflicts (pathway method + global merged)
#         8) Merged activity + bidding
#         9) Networks + permutation (3 versions, BH) + plot + save

library(NMF); library(tidyverse); library(psych); library(igraph)
library(ggraph); library(tidygraph); library(ggplot2)
library(patchwork); library(cowplot)

#======================================================================
# STEP 1: Load + filter all-zero rows 
#======================================================================
load("E:\\AID cohort\\code\\NMF\\sigLeme_widebeta.RDS")
cat("raw dim:", dim(new_Lmem_wide), "\n")
mat <- as.matrix(new_Lmem_wide)
zero_idx <- which(rowSums(mat) == 0)
cat("all-zero rows removed:", length(zero_idx), "\n")
if (length(zero_idx) > 0) {
  mat <- mat[-zero_idx, ]
  new_Lmem_wide <- new_Lmem_wide[-zero_idx, ]
}
N_total <- nrow(mat)
cat("filtered dim:", dim(mat), "N_total =", N_total, "\n")
cat("range: [", min(mat), ", ", max(mat), "]\n", sep="")

#======================================================================
# STEP 2: Robust celltype extraction
#======================================================================
ct_order <- c("DC","Mono","pDC","LDG","Baso","B","Plasma",
              "NnCD4T","NonNCD4T","NnCD8T","NonNCD8T","UTC","NK","Platelet","all")
celltype_colors <- c("DC"="#F0E442","Mono"="#F49D5C","pDC"="#B07A99","LDG"="#C95968",
                     "Baso"="#8C2522","B"="#003D81","Plasma"="#B3D1E7","NnCD4T"="#847AB3","NonNCD4T"="#C8BBD5",
                     "NnCD8T"="#F6E1EE","NonNCD8T"="#0B71AB","UTC"="#86C7B4","NK"="#96AF95","Platelet"="#B3928B","all"="#502680")
ct_sorted <- ct_order[order(nchar(ct_order), decreasing = TRUE)]
extract_celltype <- function(x, ct_sorted) {
  for (ct in ct_sorted) {
    if (grepl(paste0("Factor_\\d+[ _]", ct, "[ _]"), x)) return(ct)
  }
  return("unknown")
}
factor_meta <- data.frame(
  factor_id = rownames(mat),
  is_inhibit = grepl("^Inhibit_", rownames(mat)),
  stringsAsFactors = FALSE)
factor_meta$base_name <- sub("^Inhibit_", "", factor_meta$factor_id)
factor_meta$celltype <- sapply(factor_meta$base_name, extract_celltype, ct_sorted = ct_sorted)
cat("celltype distribution:\n"); print(table(factor_meta$celltype))
cat("unknown:", sum(factor_meta$celltype == "unknown"), "\n")
cat("inhibit:", sum(factor_meta$is_inhibit), "\n")

#======================================================================
# STEP 3: NMF rank estimation + decomposition
#======================================================================
set.seed(123456)
estim <- nmf(mat, 2:16, nrun = 50)
plot(estim)
k <- 12
nmf_result <- nmf(mat, rank = k, nrun = 50, seed = 123456)
W <- basis(nmf_result)   # 174 factor x 12 module
H <- coef(nmf_result)    # 12 module x 16 disease
colnames(W) <- paste0("M", 1:k); rownames(H) <- paste0("M", 1:k)
cat("W dim:", dim(W), "H dim:", dim(H), "\n")
#save(estim, nmf_result, file = "E:\\AID cohort\\code\\NMF\\0902_K12NMF.Rdata")

#======================================================================
# STEP 4: spearman correlation + specificity matrix
#======================================================================
cor_test_result <- psych::corr.test(t(mat), method = "spearman", adjust = "fdr")
cor_matrix <- cor_test_result$r; fdr_matrix <- cor_test_result$p
calc_specificity_matrix <- function(R) {
  s <- nrow(R); spec_mat <- matrix(NA, s, s)
  rownames(spec_mat) <- rownames(R); colnames(spec_mat) <- colnames(R)
  for (i in 1:(s-1)) for (j in (i+1):s) {
    spec <- sum(c(R[i,-j], R[-i,j]) <= R[i,j]) / (2*s - 2)
    spec_mat[i,j] <- spec; spec_mat[j,i] <- spec
  }
  diag(spec_mat) <- NA; return(spec_mat)
}
spec_matrix <- calc_specificity_matrix(cor_matrix)
cat("cor dim:", dim(cor_matrix), "fdr dim:", dim(fdr_matrix), "spec dim:", dim(spec_matrix), "\n")


#======================================================================
# STEP 5: SENSITIVITY ANALYSIS — determine top_n (BH-corrected internally)
#======================================================================
sensitivity_topn <- function(mat, W, spec_matrix, cor_matrix, fdr_matrix,
                             factor_meta, top_n_range=1:50,
                             cor_cutoff=0.3, fdr_cutoff=0.05,
                             N=nrow(cor_matrix), n_perm=200) {
  modules <- colnames(W); all_factors <- rownames(W)
  act_factors <- all_factors[!grepl("^Inhibit_",all_factors)]
  inh_factors <- all_factors[grepl("^Inhibit_",all_factors)]
  pass_corfdr <- matrix(FALSE,nrow(cor_matrix),ncol(cor_matrix))
  rownames(pass_corfdr)<-rownames(cor_matrix);colnames(pass_corfdr)<-colnames(cor_matrix)
  for(i in 1:(nrow(cor_matrix)-1)) for(j in (i+1):ncol(cor_matrix))
    if(!is.na(cor_matrix[i,j])&&!is.na(fdr_matrix[i,j]))
      if(abs(cor_matrix[i,j])>=cor_cutoff&&fdr_matrix[i,j]<=fdr_cutoff)
        pass_corfdr[i,j]<-pass_corfdr[j,i]<-TRUE
  results <- list(); idx <- 0
  for(tn in top_n_range) {
    for(type in c("activation","inhibition")) {
      fsub <- if(type=="activation") act_factors else inh_factors
      tmp <- list()
      for(mi in seq_along(modules)) {
        m <- modules[mi]
        wv <- W[fsub,m]; names(wv)<-fsub
        nt <- min(tn,length(wv)); tf <- names(sort(wv,decreasing=TRUE))[1:nt]
        n <- length(tf)
        if(n<2){tmp[[mi]]<-data.frame(top_n=tn,type=type,module=m,n_top=n,n_nodes=0,n_edges=0,conn=0,p_perm=NA,stringsAsFactors=FALSE);next}
        sc <- 1-((n-1)*(2^(-1/n)))/N
        ed <- data.frame(f1=character(),f2=character(),stringsAsFactors=FALSE)
        cb <- combn(tf,2)
        for(kk in 1:ncol(cb)){fi<-cb[1,kk];fj<-cb[2,kk]
        if(pass_corfdr[fi,fj]&&!is.na(spec_matrix[fi,fj])&&spec_matrix[fi,fj]>=sc)
          ed<-rbind(ed,data.frame(f1=fi,f2=fj,stringsAsFactors=FALSE))}
        if(nrow(ed)==0){tmp[[mi]]<-data.frame(top_n=tn,type=type,module=m,n_top=n,n_nodes=0,n_edges=0,conn=0,p_perm=NA,stringsAsFactors=FALSE);next}
        gf <- graph_from_data_frame(ed,directed=FALSE);cp <- components(gf)
        g <- induced_subgraph(gf,which(cp$membership==which.max(cp$csize)))
        n0<-gorder(g);e0<-gsize(g);mx<-n0*(n0-1)/2;cs0<-if(mx==0)0 else e0/mx
        ps <- numeric(n_perm)
        for(p in 1:n_perm){sel<-sample(all_factors,n0);sct<-1-((n0-1)*(2^(-1/n0)))/N
        se<-data.frame(f1=character(),f2=character(),stringsAsFactors=FALSE)
        if(length(sel)>=2){cs2<-combn(sel,2)
        for(kk in 1:ncol(cs2)){fi<-cs2[1,kk];fj<-cs2[2,kk]
        if(pass_corfdr[fi,fj]&&!is.na(spec_matrix[fi,fj])&&spec_matrix[fi,fj]>=sct)
          se<-rbind(se,data.frame(f1=fi,f2=fj,stringsAsFactors=FALSE))}}
        if(nrow(se)>0){gs2<-graph_from_data_frame(se,directed=FALSE);cp2<-components(gs2)
        gs3<-induced_subgraph(gs2,which(cp2$membership==which.max(cp2$csize)))
        ps[p]<-gsize(gs3)/mx}else{ps[p]<-0}}
        pp<-sum(ps>=cs0)/n_perm
        tmp[[mi]]<-data.frame(top_n=tn,type=type,module=m,n_top=n,n_nodes=n0,n_edges=e0,conn=cs0,p_perm=pp,stringsAsFactors=FALSE)
      }
      tmp_df <- do.call(rbind,tmp)
      tmp_df$p_bh <- p.adjust(tmp_df$p_perm,method="BH")
      for(mi in seq_along(modules)){idx<-idx+1;results[[idx]]<-tmp_df[mi,]}
    }
    cat(sprintf("top_n=%d done\n",tn))
  }
  return(do.call(rbind,results))
}

sens_res <- sensitivity_topn(
  mat=mat, W=W, spec_matrix=spec_matrix,
  cor_matrix=cor_matrix, fdr_matrix=fdr_matrix,
  factor_meta=factor_meta,
  top_n_range=seq(5,50,by=5),
  cor_cutoff=0.3, fdr_cutoff=0.05,
  N=nrow(cor_matrix), n_perm=10000)

sens_summary <- sens_res %>%
  group_by(top_n, type) %>%
  summarise(mean_conn=mean(conn,na.rm=TRUE),
            n_sig=sum(p_bh<0.05,na.rm=TRUE),
            .groups="drop")

plot_sensitivity <- function(sens_summary, colors=c("activation"="#E64B35","inhibition"="#4DBBD5")) {
  th <- theme_bw()+theme(
    panel.grid.major=element_line(linetype="dashed",color="grey85"),
    panel.grid.minor=element_line(linetype="dashed",color="grey92"),
    axis.text=element_text(color="black"),axis.ticks=element_line(color="black"))
  p1 <- ggplot(sens_summary,aes(top_n,mean_conn,color=type,group=type))+
    geom_line(linewidth=1)+geom_point(size=2)+
    geom_vline(data=sens_summary%>%group_by(type)%>%filter(mean_conn==max(mean_conn)),
               aes(xintercept=top_n,color=type),linetype="dashed",alpha=0.5,show.legend=FALSE)+
    scale_color_manual(values=colors)+th+
    labs(x="top_n",y="mean connection score",color=NULL)
  p2 <- ggplot(sens_summary,aes(top_n,n_sig,color=type,group=type))+
    geom_line(linewidth=1)+geom_point(size=2)+
    geom_hline(yintercept=12,linetype="dashed",color="gray")+
    scale_color_manual(values=colors)+th+
    labs(x="top_n",y="modules with p<0.05",color=NULL)
  p1/p2+plot_layout(ncol=1,guides="collect")
}

p <- plot_sensitivity(sens_summary, colors=c("activation"="#B07A99","inhibition"="#0B71AB"))
print(p)

# determine peak top_n for each type
sens_summary[sens_summary$type=="activation",]%>%filter(mean_conn==max(mean_conn))
sens_summary[sens_summary$type=="inhibition",]%>%filter(mean_conn==max(mean_conn))

#======================================================================
# STEP 6: Build module edge lists — SEPARATE top_n for act/inh
#======================================================================
build_module_edge_lists_sep <- function(mat, W, spec_matrix, cor_matrix, fdr_matrix,
                                        factor_meta = NULL,
                                        top_n_act = 25, top_n_inh = 25,
                                        cor_cutoff = 0.3, fdr_cutoff = 0.05) {
  modules <- colnames(W)
  all_factors <- rownames(W)
  act_factors <- all_factors[!grepl("^Inhibit_", all_factors)]
  inh_factors <- all_factors[grepl("^Inhibit_", all_factors)]
  cat(sprintf("activation: %d factors, top_n=%d | inhibition: %d factors, top_n=%d\n",
              length(act_factors), top_n_act, length(inh_factors), top_n_inh))
  
  build_one <- function(factor_subset, top_n, label) {
    edge_lists <- list()
    for (m in modules) {
      w_vec <- W[factor_subset, m]; names(w_vec) <- factor_subset
      n_take <- min(top_n, length(w_vec))
      top_factors <- names(sort(w_vec, decreasing = TRUE))[1:n_take]
      n <- length(top_factors)
      spec_cutoff <- 1 - ((n - 1) * (2 ^ (-1/n))) / N_total
      
      if (!is.null(factor_meta)) {
        node_ct <- setNames(factor_meta$celltype[match(top_factors, factor_meta$factor_id)], top_factors)
      } else {
        node_ct <- setNames(rep(NA, length(top_factors)), top_factors)
      }
      
      edges <- data.frame(factor1=character(), factor2=character(),
                          factor1_celltype=character(), factor2_celltype=character(),
                          correlation=numeric(), fdr=numeric(), specificity=numeric(),
                          spec_cutoff=numeric(), stringsAsFactors=FALSE)
      if (length(top_factors) >= 2) {
        combs <- combn(top_factors, 2)
        for (kk in 1:ncol(combs)) {
          fi <- combs[1,kk]; fj <- combs[2,kk]
          rij <- cor_matrix[fi,fj]; fdr_ij <- fdr_matrix[fi,fj]; spec_ij <- spec_matrix[fi,fj]
          if (!is.na(rij) && !is.na(fdr_ij) && !is.na(spec_ij))
            if (abs(rij) >= cor_cutoff && fdr_ij <= fdr_cutoff && spec_ij >= spec_cutoff)
              edges <- rbind(edges, data.frame(factor1=fi, factor2=fj,
                                               factor1_celltype=node_ct[fi], factor2_celltype=node_ct[fj],
                                               correlation=rij, fdr=fdr_ij, specificity=spec_ij,
                                               spec_cutoff=spec_cutoff, stringsAsFactors=FALSE))
        }
      }
      cat(sprintf("---%s [%s]: top %d factors, %d edges (spec_cut=%.3f)---\n",
                  m, label, n, nrow(edges), spec_cutoff))
      edge_lists[[m]] <- edges
    }
    return(edge_lists)
  }
  
  act_edges <- build_one(act_factors, top_n_act, "activation")
  inh_edges <- build_one(inh_factors, top_n_inh, "inhibition")
  return(list(activation = act_edges, inhibition = inh_edges))
}


# peak top_n from STEP 5 sensitivity (group-specific max)
act_df <- sens_summary %>% filter(type=="activation")
inh_df <- sens_summary %>% filter(type=="inhibition")
peak_act <- act_df$top_n[which.max(act_df$mean_conn)]
peak_inh <- inh_df$top_n[which.max(inh_df$mean_conn)]
cat(sprintf("Using peak top_n: activation=%d, inhibition=%d\n", peak_act, peak_inh))

edges_sep <- build_module_edge_lists_sep(
  mat = mat, W = W, spec_matrix = spec_matrix,
  cor_matrix = cor_matrix, fdr_matrix = fdr_matrix,
  factor_meta = factor_meta,
  top_n_act = peak_act,
  top_n_inh = peak_inh,
  cor_cutoff = 0.3, fdr_cutoff = 0.05)

module_edges_act <- edges_sep$activation
module_edges_inh <- edges_sep$inhibition
cat("\n=== activation edges per module ===\n"); print(lapply(module_edges_act, nrow))
cat("\n=== inhibition edges per module ===\n"); print(lapply(module_edges_inh, nrow))

#======================================================================
# STEP 7: Module activity + disease assignment — SEPARATE act/inh
#   Two assignment methods: (1) which.max, (2) bidding
#======================================================================
module_activity <- function(mat, W, module_edges) {
  diseases <- colnames(mat); modules <- names(module_edges)
  act_mat <- matrix(0, length(modules), length(diseases))
  rownames(act_mat) <- modules; colnames(act_mat) <- diseases
  for (j in seq_along(modules)) {
    edges <- module_edges[[modules[j]]]
    if (nrow(edges) == 0) next
    nodes <- intersect(unique(c(edges$factor1, edges$factor2)), rownames(mat))
    if (length(nodes) == 0) next
    for (i in seq_along(diseases))
      act_mat[j, i] <- sum(mat[nodes, diseases[i]] * W[nodes, modules[j]])
  }
  q99 <- as.numeric(quantile(as.vector(act_mat), 0.99))
  scaled <- act_mat
  if (q99 > 0) { scaled[act_mat > q99] <- 1; scaled[act_mat <= q99] <- act_mat[act_mat <= q99] / q99 }
  mod_class <- data.frame(disease = diseases,
                          module = apply(scaled, 2, function(x) rownames(scaled)[which.max(x)]),
                          stringsAsFactors = FALSE)
  return(list(activity = as.data.frame(scaled), class = mod_class))
}

# --- method 1: which.max assignment ---
mod_act <- module_activity(mat = mat, W = W, module_edges = module_edges_act)
mod_inh <- module_activity(mat = mat, W = W, module_edges = module_edges_inh)
cat("\n=== activation (which.max) ===\n"); print(mod_act$class)
cat("\n=== inhibition (which.max) ===\n"); print(mod_inh$class)

# --- method 2: bidding assignment ---
assign_modules_bidding <- function(activity_mat, share_threshold = 1.0) {
  diseases <- colnames(activity_mat)
  modules  <- rownames(activity_mat)
  prefs <- lapply(diseases, function(d) names(sort(activity_mat[, d], decreasing = TRUE)))
  names(prefs) <- diseases
  mod_assigned <- setNames(vector("list", length(modules)), modules)
  disease_done <- setNames(rep(FALSE, length(diseases)), diseases)
  for (rnd in 1:(length(diseases) + 2)) {
    progress <- FALSE
    for (d in diseases[!disease_done]) {
      for (m in prefs[[d]]) {
        d_act   <- activity_mat[m, d]
        current <- mod_assigned[[m]]
        if (length(current) == 0) {
          mod_assigned[[m]] <- d; disease_done[d] <- TRUE; progress <- TRUE; break
        }
        occupant_acts <- sapply(current, function(od) activity_mat[m, od])
        if (d_act >= share_threshold && all(occupant_acts >= share_threshold)) {
          mod_assigned[[m]] <- c(current, d); disease_done[d] <- TRUE; progress <- TRUE; break
        }
        if (d_act > max(occupant_acts)) {
          for (od in current) disease_done[od] <- FALSE
          mod_assigned[[m]] <- d; disease_done[d] <- TRUE; progress <- TRUE; break
        }
      }
    }
    if (!progress || all(disease_done)) break
  }
  n_kicked <- sum(!disease_done)
  for (d in diseases[!disease_done]) {
    best_m <- names(sort(activity_mat[, d], decreasing = TRUE))[1]
    mod_assigned[[best_m]] <- c(mod_assigned[[best_m]], d)
    disease_done[d] <- TRUE
  }
  res <- data.frame(disease = diseases, module = NA_character_,
                    activity = NA_real_, stringsAsFactors = FALSE)
  for (m in modules) {
    for (d in mod_assigned[[m]]) {
      res$module[res$disease == d]   <- m
      res$activity[res$disease == d] <- activity_mat[m, d]
    }
  }
  res <- res[order(res$module, -res$activity), ]; rownames(res) <- NULL
  cat(sprintf("[bidding] %d/%d modules used, %d diseases assigned\n",
              sum(sapply(mod_assigned, length) > 0), length(modules), nrow(res)))
  return(list(module_assigned = mod_assigned, result = res))
}

bidding_act <- assign_modules_bidding(activity_mat = as.matrix(mod_act$activity))
bidding_inh <- assign_modules_bidding(activity_mat = as.matrix(mod_inh$activity))
bidding_act$result <- bidding_act$result[bidding_act$result$activity != 0, ]
bidding_inh$result <- bidding_inh$result[bidding_inh$result$activity != 0, ]
cat("\n=== activation (bidding) ===\n"); print(bidding_act$result)
cat("\n=== inhibition (bidding) ===\n"); print(bidding_inh$result)

#======================================================================
# STEP 8: Resolve conflicting factors — TWO methods
#   Method 1 (pathway): per-pathway beta-mean, returns act+inh separate
#   Method 2 (global):  merged connected-component, returns merged edges
#======================================================================

#---- Method 1: per-pathway beta-mean resolution ----
resolve_conflicts_pathway <- function(bidding_act, bidding_inh,
                                      module_edges_act, module_edges_inh,
                                      mat, factor_meta, n_core_words = 3) {
  get_core <- function(fname) {
    base <- sub("^(Inhibit_)?Factor_\\d+_", "", fname)
    ct <- factor_meta$celltype[match(fname, factor_meta$factor_id)]
    pw <- sub(paste0("^", ct, "_"), "", base)
    paste(head(strsplit(pw, "_")[[1]], n_core_words), collapse="_")
  }
  diseases <- intersect(bidding_act$result$disease, bidding_inh$result$disease)
  rem_am <- setNames(vector("list", length(module_edges_act)), names(module_edges_act))
  rem_im <- setNames(vector("list", length(module_edges_inh)), names(module_edges_inh))
  for (d in diseases) {
    am <- bidding_act$result$module[bidding_act$result$disease == d]
    im <- bidding_inh$result$module[bidding_inh$result$disease == d]
    ae <- module_edges_act[[am]]; ie <- module_edges_inh[[im]]
    an <- unique(c(ae$factor1, ae$factor2)); in_ <- unique(c(ie$factor1, ie$factor2))
    ac <- sapply(an, get_core); ic <- sapply(in_, get_core)
    ra <- character(); ri <- character()
    for (core in intersect(unique(ac), unique(ic))) {
      af <- names(ac[ac == core]); if_ <- names(ic[ic == core])
      are <- ae[ae$factor1 %in% af | ae$factor2 %in% af, ]
      ire <- ie[ie$factor1 %in% if_ | ie$factor2 %in% if_, ]
      arf <- unique(c(are$factor1, are$factor2)); irf <- unique(c(ire$factor1, ire$factor2))
      if (abs(mean(mat[arf, d], na.rm=TRUE)) >= abs(mean(mat[irf, d], na.rm=TRUE)))
        ri <- c(ri, irf) else ra <- c(ra, arf)
    }
    rem_am[[am]] <- union(rem_am[[am]], ra); rem_im[[im]] <- union(rem_im[[im]], ri)
  }
  fa <- lapply(names(module_edges_act), function(m) {
    e <- module_edges_act[[m]]; r <- rem_am[[m]]
    e[!e$factor1 %in% r & !e$factor2 %in% r, ]
  }); names(fa) <- names(module_edges_act)
  fi <- lapply(names(module_edges_inh), function(m) {
    e <- module_edges_inh[[m]]; r <- rem_im[[m]]
    e[!e$factor1 %in% r & !e$factor2 %in% r, ]
  }); names(fi) <- names(module_edges_inh)
  return(list(activation = fa, inhibition = fi))
}

#---- Method 2: merged connected-component resolution ----
resolve_conflicts_global <- function(bidding_act, bidding_inh,
                                     module_edges_act, module_edges_inh,
                                     mat, factor_meta) {
  diseases <- intersect(bidding_act$result$disease, bidding_inh$result$disease)
  rem_am <- setNames(vector("list", length(module_edges_act)), names(module_edges_act))
  rem_im <- setNames(vector("list", length(module_edges_inh)), names(module_edges_inh))
  for (d in diseases) {
    am <- bidding_act$result$module[bidding_act$result$disease == d]
    im <- bidding_inh$result$module[bidding_inh$result$disease == d]
    ae <- module_edges_act[[am]]; ie <- module_edges_inh[[im]]
    af <- unique(c(ae$factor1, ae$factor2)); if_ <- unique(c(ie$factor1, ie$factor2))
    g <- graph_from_data_frame(rbind(ae[,1:2], ie[,1:2]), directed=FALSE)
    comp <- components(g)
    ra <- character(); ri <- character()
    for (ci in 1:comp$no) {
      cn <- names(comp$membership[comp$membership == ci])
      ca <- intersect(cn, af); ci_ <- intersect(cn, if_)
      if (length(ca) == 0 || length(ci_) == 0) next
      if (abs(mean(mat[ca, d], na.rm=TRUE)) >= abs(mean(mat[ci_, d], na.rm=TRUE)))
        ri <- c(ri, ci_) else ra <- c(ra, ca)
    }
    rem_am[[am]] <- union(rem_am[[am]], ra); rem_im[[im]] <- union(rem_im[[im]], ri)
  }
  fa <- lapply(names(module_edges_act), function(m) {
    e <- module_edges_act[[m]]; r <- rem_am[[m]]
    e[!e$factor1 %in% r & !e$factor2 %in% r, ]
  }); names(fa) <- names(module_edges_act)
  fi <- lapply(names(module_edges_inh), function(m) {
    e <- module_edges_inh[[m]]; r <- rem_im[[m]]
    e[!e$factor1 %in% r & !e$factor2 %in% r, ]
  }); names(fi) <- names(module_edges_inh)
  merged <- lapply(names(fa), function(m) rbind(fa[[m]], fi[[m]]))
  names(merged) <- names(fa)
  return(merged)
}

res_pw <- resolve_conflicts_pathway(
  bidding_act = bidding_act, bidding_inh = bidding_inh,
  module_edges_act = module_edges_act, module_edges_inh = module_edges_inh,
  mat = mat, factor_meta = factor_meta, n_core_words = 3)

module_edges_act_pw  <- res_pw$activation
module_edges_inh_pw  <- res_pw$inhibition
module_edges_merged  <- resolve_conflicts_global(
  bidding_act = bidding_act, bidding_inh = bidding_inh,
  module_edges_act = module_edges_act, module_edges_inh = module_edges_inh,
  mat = mat, factor_meta = factor_meta)

# check edges per module
lapply(module_edges_act_pw, nrow)
lapply(module_edges_inh_pw, nrow)
lapply(module_edges_merged, nrow)


mod_merged <- module_activity(mat = mat, W = W, module_edges = module_edges_merged)
bidding_merged <- assign_modules_bidding(activity_mat = as.matrix(mod_merged$activity))
bidding_merged$result <- bidding_merged$result[bidding_merged$result$activity != 0, ]
cat("\n=== merged (bidding) ===\n"); print(bidding_merged$result)#split module for each disease again



#======================================================================
# STEP 9: Network score + permutation test function  
#======================================================================
gray_gradient <- colorRampPalette(c("#EBEBEB", "#D7D7D7", "#B4B4B4", "#8E8E8E", "#5F5F5F"))(100)

network_score_p <- function(module_edges, cor_matrix, fdr_matrix, spec_matrix,
                            factor_meta, celltype_colors,
                            cor_cutoff = 0.3, fdr_cutoff = 0.05,
                            N = nrow(cor_matrix), n_perm = 1000,
                            node_size = 7) {
  network_list <- list()
  pass_corfdr <- matrix(FALSE, nrow(cor_matrix), ncol(cor_matrix))
  rownames(pass_corfdr) <- rownames(cor_matrix); colnames(pass_corfdr) <- colnames(cor_matrix)
  for (i in 1:(nrow(cor_matrix)-1)) {
    for (j in (i+1):ncol(cor_matrix)) {
      if (!is.na(cor_matrix[i,j]) && !is.na(fdr_matrix[i,j]))
        if (abs(cor_matrix[i,j]) >= cor_cutoff && fdr_matrix[i,j] <= fdr_cutoff)
          pass_corfdr[i,j] <- pass_corfdr[j,i] <- TRUE
    }
  }
  
  for (m in names(module_edges)) {
    edges <- module_edges[[m]]
    if (nrow(edges) == 0) { cat(sprintf("%s: no edges, skip\n", m)); next }
    
    g_full <- graph_from_data_frame(edges[, c("factor1","factor2")], directed = FALSE)
    comp <- components(g_full)
    g <- induced_subgraph(g_full, which(comp$membership == which.max(comp$csize)))
    
    edge_ends <- ends(g, E(g))
    E(g)$specificity <- apply(edge_ends, 1, function(x) {
      idx <- which((edges$factor1 == x[1] & edges$factor2 == x[2]) |
                     (edges$factor1 == x[2] & edges$factor2 == x[1]))
      if (length(idx) == 1) return(edges$specificity[idx]) else return(NA)
    })
    
    node_names <- V(g)$name
    node_ct <- factor_meta$celltype[match(node_names, factor_meta$factor_id)]
    node_fill <- celltype_colors[match(node_ct, names(celltype_colors))]
    node_fill[is.na(node_fill)] <- "#888888"
    V(g)$fill <- node_fill
    V(g)$shape <- ifelse(grepl("^Inhibit_", node_names), 16, 21)
    V(g)$color <- ifelse(grepl("^Inhibit_", node_names), node_fill, "black")
    V(g)$label <- ifelse(grepl("^Inhibit_", node_names),
                         paste0("iF", sub("^Inhibit_Factor_(\\d+).*", "\\1", node_names)),
                         paste0("F", sub("^Factor_(\\d+).*", "\\1", node_names)))
    
    n_nodes0 <- gorder(g); n_edges0 <- gsize(g)
    max_edges0 <- n_nodes0 * (n_nodes0 - 1) / 2
    conn_score0 <- if (max_edges0 == 0) 0 else n_edges0 / max_edges0
    
    perm_scores <- numeric(n_perm)
    all_factors <- rownames(cor_matrix)
    for (p in 1:n_perm) {
      sel <- sample(all_factors, n_nodes0)
      spec_cut <- 1 - ((n_nodes0 - 1) * (2 ^ (-1/n_nodes0))) / N
      sel_edges <- data.frame(f1=character(), f2=character(), stringsAsFactors=FALSE)
      if (length(sel) >= 2) {
        combs <- combn(sel, 2)
        for (kk in 1:ncol(combs)) {
          fi <- combs[1,kk]; fj <- combs[2,kk]
          if (pass_corfdr[fi,fj] && !is.na(spec_matrix[fi,fj]) && spec_matrix[fi,fj] >= spec_cut)
            sel_edges <- rbind(sel_edges, data.frame(f1=fi, f2=fj, stringsAsFactors=FALSE))
        }
      }
      if (nrow(sel_edges) > 0) {
        g_sel <- graph_from_data_frame(sel_edges, directed = FALSE)
        comp_sel <- components(g_sel)
        g_sel <- induced_subgraph(g_sel, which(comp_sel$membership == which.max(comp_sel$csize)))
        perm_scores[p] <- gsize(g_sel) / max_edges0
      } else { perm_scores[p] <- 0 }
    }
    p_perm <- sum(perm_scores >= conn_score0) / n_perm
    
    tg <- as_tbl_graph(g)
    sig_label <- ifelse(p_perm < 0.001, "***",
                        ifelse(p_perm < 0.01, "**",
                               ifelse(p_perm < 0.05, "*", "ns")))
    circle_plot <- ggraph(tg, layout = "circle") +
      geom_edge_link(aes(color = specificity), alpha = 0.7, show.legend = TRUE) +
      scale_edge_color_gradientn(colors = gray_gradient, name = "Specificity") +
      geom_node_point(size = node_size, stroke = 0.8,
                      aes(fill = fill, color = color, shape = shape),
                      show.legend = FALSE) +
      scale_fill_identity() + scale_color_identity() + scale_shape_identity() +
      geom_node_text(aes(label = label), size = 3, nudge_x = 0.12, nudge_y = -0.12) +
      theme_void() +
      ggtitle(sprintf("Module %s (nodes=%d, edges=%d, conn=%.3f, p=%.4f %s)",
                      sub("M","",m), n_nodes0, n_edges0, conn_score0, p_perm, sig_label)) +
      theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 10))
    
    edge_info <- data.frame(factor1 = edge_ends[,1], factor2 = edge_ends[,2],
                            correlation = NA, fdr = NA, specificity = E(g)$specificity,
                            stringsAsFactors = FALSE)
    for (i in 1:nrow(edge_info)) {
      idx <- which((edges$factor1 == edge_info$factor1[i] & edges$factor2 == edge_info$factor2[i]) |
                     (edges$factor1 == edge_info$factor2[i] & edges$factor2 == edge_info$factor1[i]))
      if (length(idx) == 1) {
        edge_info$correlation[i] <- edges$correlation[idx]
        edge_info$fdr[i] <- edges$fdr[idx]
      }
    }
    
    cat(sprintf("%s: nodes=%d, edges=%d, conn=%.3f, p=%.4f %s\n",
                m, n_nodes0, n_edges0, conn_score0, p_perm, sig_label))
    network_list[[m]] <- list(igraph = circle_plot, graph = g,
                              connection_score = conn_score0, permutation_p = p_perm,
                              n_nodes = n_nodes0, n_edges = n_edges0,
                              edges_in_graph = edge_info)
  }
  return(network_list)
}
# fill and width of line dependent of specificity

network_score_p <- function(module_edges, cor_matrix, fdr_matrix, spec_matrix,
                            factor_meta, celltype_colors,
                            cor_cutoff = 0.3, fdr_cutoff = 0.05,
                            N = nrow(cor_matrix), n_perm = 1000,
                            node_size = 9, edge_width = 1) {
  network_list <- list()
  pass_corfdr <- matrix(FALSE, nrow(cor_matrix), ncol(cor_matrix))
  rownames(pass_corfdr) <- rownames(cor_matrix); colnames(pass_corfdr) <- colnames(cor_matrix)
  for (i in 1:(nrow(cor_matrix)-1)) {
    for (j in (i+1):ncol(cor_matrix)) {
      if (!is.na(cor_matrix[i,j]) && !is.na(fdr_matrix[i,j]))
        if (abs(cor_matrix[i,j]) >= cor_cutoff && fdr_matrix[i,j] <= fdr_cutoff)
          pass_corfdr[i,j] <- pass_corfdr[j,i] <- TRUE
    }
  }
  
  for (m in names(module_edges)) {
    edges <- module_edges[[m]]
    if (nrow(edges) == 0) { cat(sprintf("%s: no edges, skip\n", m)); next }
    
    g_full <- graph_from_data_frame(edges[, c("factor1","factor2")], directed = FALSE)
    comp <- components(g_full)
    g <- induced_subgraph(g_full, which(comp$membership == which.max(comp$csize)))
    
    edge_ends <- ends(g, E(g))
    E(g)$specificity <- apply(edge_ends, 1, function(x) {
      idx <- which((edges$factor1 == x[1] & edges$factor2 == x[2]) |
                     (edges$factor1 == x[2] & edges$factor2 == x[1]))
      if (length(idx) == 1) return(edges$specificity[idx]) else return(NA)
    })
    
    node_names <- V(g)$name
    node_ct <- factor_meta$celltype[match(node_names, factor_meta$factor_id)]
    node_fill <- celltype_colors[match(node_ct, names(celltype_colors))]
    node_fill[is.na(node_fill)] <- "#888888"
    V(g)$fill <- node_fill
    V(g)$shape <- ifelse(grepl("^Inhibit_", node_names), 16, 21)
    V(g)$color <- ifelse(grepl("^Inhibit_", node_names), node_fill, "black")
    V(g)$label <- ifelse(grepl("^Inhibit_", node_names),
                         paste0("iF", sub("^Inhibit_Factor_(\\d+).*", "\\1", node_names)),
                         paste0("F", sub("^Factor_(\\d+).*", "\\1", node_names)))
    
    n_nodes0 <- gorder(g); n_edges0 <- gsize(g)
    max_edges0 <- n_nodes0 * (n_nodes0 - 1) / 2
    conn_score0 <- if (max_edges0 == 0) 0 else n_edges0 / max_edges0
    
    perm_scores <- numeric(n_perm)
    all_factors <- rownames(cor_matrix)
    for (p in 1:n_perm) {
      sel <- sample(all_factors, n_nodes0)
      spec_cut <- 1 - ((n_nodes0 - 1) * (2 ^ (-1/n_nodes0))) / N
      sel_edges <- data.frame(f1=character(), f2=character(), stringsAsFactors=FALSE)
      if (length(sel) >= 2) {
        combs <- combn(sel, 2)
        for (kk in 1:ncol(combs)) {
          fi <- combs[1,kk]; fj <- combs[2,kk]
          if (pass_corfdr[fi,fj] && !is.na(spec_matrix[fi,fj]) && spec_matrix[fi,fj] >= spec_cut)
            sel_edges <- rbind(sel_edges, data.frame(f1=fi, f2=fj, stringsAsFactors=FALSE))
        }
      }
      if (nrow(sel_edges) > 0) {
        g_sel <- graph_from_data_frame(sel_edges, directed = FALSE)
        comp_sel <- components(g_sel)
        g_sel <- induced_subgraph(g_sel, which(comp_sel$membership == which.max(comp_sel$csize)))
        perm_scores[p] <- gsize(g_sel) / max_edges0
      } else { perm_scores[p] <- 0 }
    }
    p_perm <- sum(perm_scores >= conn_score0) / n_perm
    
    tg <- as_tbl_graph(g)
    sig_label <- ifelse(p_perm < 0.001, "***",
                        ifelse(p_perm < 0.01, "**",
                               ifelse(p_perm < 0.05, "*", "ns")))
    circle_plot <- ggraph(tg, layout = "circle") +
      geom_edge_link(aes(color = specificity), alpha = 0.7,
                     width = edge_width, show.legend = TRUE) +
      scale_edge_color_gradientn(colors = gray_gradient, name = "Specificity") +
      geom_node_point(size = node_size, stroke = 0.8,
                      aes(fill = fill, color = color, shape = shape),
                      show.legend = FALSE) +
      scale_fill_identity() + scale_color_identity() + scale_shape_identity() +
      geom_node_text(aes(label = label), size = 3, nudge_x = 0.12, nudge_y = -0.12) +
      theme_void() +
      ggtitle(sprintf("Module %s (nodes=%d, edges=%d, conn=%.3f, p=%.4f %s)",
                      sub("M","",m), n_nodes0, n_edges0, conn_score0, p_perm, sig_label)) +
      theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 10))
    
    edge_info <- data.frame(factor1 = edge_ends[,1], factor2 = edge_ends[,2],
                            correlation = NA, fdr = NA, specificity = E(g)$specificity,
                            stringsAsFactors = FALSE)
    for (i in 1:nrow(edge_info)) {
      idx <- which((edges$factor1 == edge_info$factor1[i] & edges$factor2 == edge_info$factor2[i]) |
                     (edges$factor1 == edge_info$factor2[i] & edges$factor2 == edge_info$factor1[i]))
      if (length(idx) == 1) {
        edge_info$correlation[i] <- edges$correlation[idx]
        edge_info$fdr[i] <- edges$fdr[idx]
      }
    }
    
    cat(sprintf("%s: nodes=%d, edges=%d, conn=%.3f, p=%.4f %s\n",
                m, n_nodes0, n_edges0, conn_score0, p_perm, sig_label))
    network_list[[m]] <- list(igraph = circle_plot, graph = g,
                              connection_score = conn_score0, permutation_p = p_perm,
                              n_nodes = n_nodes0, n_edges = n_edges0,
                              edges_in_graph = edge_info)
  }
  return(network_list)
}
#width of line dependent of specificity

network_score_p <- function(module_edges, cor_matrix, fdr_matrix, spec_matrix,
                            factor_meta, celltype_colors,
                            cor_cutoff = 0.3, fdr_cutoff = 0.05,
                            N = nrow(cor_matrix), n_perm = 1000,
                            nodetext_size=3.5,
                            node_size = 9, edge_width = 0.5) {
  network_list <- list()
  pass_corfdr <- matrix(FALSE, nrow(cor_matrix), ncol(cor_matrix))
  rownames(pass_corfdr) <- rownames(cor_matrix); colnames(pass_corfdr) <- colnames(cor_matrix)
  for (i in 1:(nrow(cor_matrix)-1)) {
    for (j in (i+1):ncol(cor_matrix)) {
      if (!is.na(cor_matrix[i,j]) && !is.na(fdr_matrix[i,j]))
        if (abs(cor_matrix[i,j]) >= cor_cutoff && fdr_matrix[i,j] <= fdr_cutoff)
          pass_corfdr[i,j] <- pass_corfdr[j,i] <- TRUE
    }
  }
  
  # helper: luminance of hex color (>0.5 light -> black text; <0.5 dark -> white text)
  col_lum <- function(hex) {
    rgb <- col2rgb(hex) / 255
    0.299 * rgb[1] + 0.587 * rgb[2] + 0.114 * rgb[3]
  }
  
  for (m in names(module_edges)) {
    edges <- module_edges[[m]]
    if (nrow(edges) == 0) { cat(sprintf("%s: no edges, skip\n", m)); next }
    
    g_full <- graph_from_data_frame(edges[, c("factor1","factor2")], directed = FALSE)
    comp <- components(g_full)
    g <- induced_subgraph(g_full, which(comp$membership == which.max(comp$csize)))
    
    edge_ends <- ends(g, E(g))
    E(g)$specificity <- apply(edge_ends, 1, function(x) {
      idx <- which((edges$factor1 == x[1] & edges$factor2 == x[2]) |
                     (edges$factor1 == x[2] & edges$factor2 == x[1]))
      if (length(idx) == 1) return(edges$specificity[idx]) else return(NA)
    })
    
    node_names <- V(g)$name
    node_ct <- factor_meta$celltype[match(node_names, factor_meta$factor_id)]
    node_fill <- celltype_colors[match(node_ct, names(celltype_colors))]
    node_fill[is.na(node_fill)] <- "#888888"
    V(g)$fill <- node_fill
    V(g)$shape <- ifelse(grepl("^Inhibit_", node_names), 16, 21)
    V(g)$color <- ifelse(grepl("^Inhibit_", node_names), node_fill, "black")
    V(g)$label <- ifelse(grepl("^Inhibit_", node_names),
                         paste0("iF", sub("^Inhibit_Factor_(\\d+).*", "\\1", node_names)),
                         paste0("F", sub("^Factor_(\\d+).*", "\\1", node_names)))
    # text color inside node: black on light fill, white on dark fill
    node_text_color <- ifelse(sapply(node_fill, col_lum) > 0.5, "black", "white")
    
    n_nodes0 <- gorder(g); n_edges0 <- gsize(g)
    max_edges0 <- n_nodes0 * (n_nodes0 - 1) / 2
    conn_score0 <- if (max_edges0 == 0) 0 else n_edges0 / max_edges0
    
    perm_scores <- numeric(n_perm)
    all_factors <- rownames(cor_matrix)
    for (p in 1:n_perm) {
      sel <- sample(all_factors, n_nodes0)
      spec_cut <- 1 - ((n_nodes0 - 1) * (2 ^ (-1/n_nodes0))) / N
      sel_edges <- data.frame(f1=character(), f2=character(), stringsAsFactors=FALSE)
      if (length(sel) >= 2) {
        combs <- combn(sel, 2)
        for (kk in 1:ncol(combs)) {
          fi <- combs[1,kk]; fj <- combs[2,kk]
          if (pass_corfdr[fi,fj] && !is.na(spec_matrix[fi,fj]) && spec_matrix[fi,fj] >= spec_cut)
            sel_edges <- rbind(sel_edges, data.frame(f1=fi, f2=fj, stringsAsFactors=FALSE))
        }
      }
      if (nrow(sel_edges) > 0) {
        g_sel <- graph_from_data_frame(sel_edges, directed = FALSE)
        comp_sel <- components(g_sel)
        g_sel <- induced_subgraph(g_sel, which(comp_sel$membership == which.max(comp_sel$csize)))
        perm_scores[p] <- gsize(g_sel) / max_edges0
      } else { perm_scores[p] <- 0 }
    }
    p_perm <- sum(perm_scores >= conn_score0) / n_perm
    
    tg <- as_tbl_graph(g)
    sig_label <- ifelse(p_perm < 0.001, "***",
                        ifelse(p_perm < 0.01, "**",
                               ifelse(p_perm < 0.05, "*", "ns")))
    circle_plot <- ggraph(tg, layout = "circle") +
      geom_edge_link(aes(color = specificity), alpha = 0.7,
                     width = edge_width, show.legend = TRUE) +
      scale_edge_color_gradientn(colors = gray_gradient, name = "Specificity") +
      geom_node_point(size = node_size, stroke = 0.8,
                      aes(fill = fill, color = color, shape = shape),
                      show.legend = FALSE) +
      scale_fill_identity() + scale_color_identity() + scale_shape_identity() +
      geom_node_text(aes(label = label), size = nodetext_size, color = node_text_color) +
      theme_void() +
      ggtitle(sprintf("Module %s (conn=%.3f, p.perm=%.4f %s)",
                      sub("M","",m), conn_score0, p_perm, sig_label))+
      theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 10))
    
    edge_info <- data.frame(factor1 = edge_ends[,1], factor2 = edge_ends[,2],
                            correlation = NA, fdr = NA, specificity = E(g)$specificity,
                            stringsAsFactors = FALSE)
    for (i in 1:nrow(edge_info)) {
      idx <- which((edges$factor1 == edge_info$factor1[i] & edges$factor2 == edge_info$factor2[i]) |
                     (edges$factor1 == edge_info$factor2[i] & edges$factor2 == edge_info$factor1[i]))
      if (length(idx) == 1) {
        edge_info$correlation[i] <- edges$correlation[idx]
        edge_info$fdr[i] <- edges$fdr[idx]
      }
    }
    
    cat(sprintf("%s: nodes=%d, edges=%d, conn=%.3f, p=%.4f %s\n",
                m, n_nodes0, n_edges0, conn_score0, p_perm, sig_label))
    network_list[[m]] <- list(igraph = circle_plot, graph = g,
                              connection_score = conn_score0, permutation_p = p_perm,
                              n_nodes = n_nodes0, n_edges = n_edges0,
                              edges_in_graph = edge_info)
  }
  return(network_list)
}
# factor text within point,deep color, white text, light color, black text
network_act_pw <- network_score_p(
  module_edges = module_edges_act_pw,
  cor_matrix = cor_matrix, fdr_matrix = fdr_matrix, spec_matrix = spec_matrix,
  factor_meta = factor_meta, celltype_colors = celltype_colors,
  cor_cutoff = 0.3, fdr_cutoff = 0.05,
  N = nrow(cor_matrix), n_perm = 1000, node_size = 15)

network_inh_pw <- network_score_p(
  module_edges = module_edges_inh_pw,
  cor_matrix = cor_matrix, fdr_matrix = fdr_matrix, spec_matrix = spec_matrix,
  factor_meta = factor_meta, celltype_colors = celltype_colors,
  cor_cutoff = 0.3, fdr_cutoff = 0.05,
  N = nrow(cor_matrix), n_perm = 1000, node_size = 13)

network_merged <- network_score_p(
  module_edges = module_edges_merged,
  cor_matrix = cor_matrix, fdr_matrix = fdr_matrix, spec_matrix = spec_matrix,
  factor_meta = factor_meta, celltype_colors = celltype_colors,
  cor_cutoff = 0.3, fdr_cutoff = 0.05,
  N = nrow(cor_matrix), n_perm = 1000, node_size = 15)

#---- BH correction for each version ----
p_act_pw  <- p.adjust(sapply(network_act_pw, function(x) x$permutation_p), method = "BH")
p_inh_pw  <- p.adjust(sapply(network_inh_pw, function(x) x$permutation_p), method = "BH")
#p_merged  <- p.adjust(sapply(network_merged, function(x) x$permutation_p), method = "BH")
Presult=data.frame(module = names(p_act_pw), p_act = p_act_pw, p_inh = p_inh_pw, p_merged = p_merged)
write.xlsx(Presult,file="E:/AID cohort/code/NMF/0903_NMFnetwork_permutaionPvalue.xlsx")


#---- plot ----
act_net<-plot_grid(plotlist = lapply(network_act_pw, function(x) x$igraph), ncol = 6, align = "hv")
inh_net<-plot_grid(plotlist = lapply(network_inh_pw, function(x) x$igraph), ncol = 4, align = "hv")
#merge_net<-plot_grid(plotlist = lapply(network_merged, function(x) x$igraph), ncol = 4, align = "hv")
#ggsave(merge_net,file="E:\\AID cohort\\code\\NMF\\merge_module.pdf",width = 25.59,height = 17.71)
#25.59 17.71
ggsave(act_net,file="E:\\AID cohort\\code\\NMF\\activation_module.pdf",width = 18.74,height = 5.02)
ggsave(inh_net,file="E:\\AID cohort\\code\\NMF\\inhibition_module.pdf",width = 10.47,height = 5.49)
save.image(file="E:\\AID cohort\\code\\NMF\\0902_NMF.Rdata")



###plot connection score and adjusted p (heatmap)
data1<-Presult[,c(1,2,3)]
data2<-data.frame(
c_act=as.numeric(unlist(lapply(network_act_pw,function(x){x$connection_score}))),
c_int=as.numeric(unlist(lapply(network_inh_pw,function(x){x$connection_score}))))
data<-cbind(data1,data2)
data

##(1)pheatmap
library(pheatmap); library(patchwork)
rownames(data) <- data$module
p_mat <- as.matrix(data[,c("p_act","p_inh")]); colnames(p_mat) <- c("Active","Inhibitory")
p_mat[p_mat==0] <- 1e-4
sig <- matrix(ifelse(p_mat<0.001,"***",ifelse(p_mat<0.01,"**",ifelse(p_mat<0.05,"*",""))),nrow(p_mat))
pcol <- colorRampPalette(c("#FFF5F0","#FDD0A2","#FC8D59","#D7301F","#7F0000"))(100)
ph <- pheatmap(-log10(p_mat),color=pcol,cluster_rows=F,cluster_cols=F,display_numbers=sig,
               number_color="black",fontsize_number=11,border_color="white",cellwidth=45,cellheight=22,
               main="-log10(adj. P-value)",angle_col=0,fontsize=10,legend_breaks=c(1,2,3,4),
               legend_labels=c("1","2","3",">=4"),silent=T)

c_mat <- as.matrix(data[,c("c_act","c_int")]); colnames(c_mat) <- c("Active","Inhibitory")
ccol <- colorRampPalette(c("#F7FBFF","#DEEBF7","#9ECAE1","#4292C6","#08306B"))(100)
ch <- pheatmap(c_mat,color=ccol,cluster_rows=F,cluster_cols=F,display_numbers=T,number_format="%.2f",
               number_color="black",fontsize_number=9,border_color="white",cellwidth=45,cellheight=22,
               main="Connection score",angle_col=0,fontsize=10,silent=T)

wrap_plots(ph[[4]],ch[[4]],ncol=2)

##(2) disease- module heatmap
library(tidyverse); library(ggplot2); library(patchwork)
disease_colors <- c(
  "KD" = "#2C73B9", "JIA" = "#7C7EB8", "SLE" = "#9284B7",
  "Uncategorized" = "#A9799D", "Protein homeostasis" = "#D69A9D",
  "TBK1_IRF3" = "#8F7966", "Cytoskeleton and small GTPase" = "#A6A6A6",
  "Negative regulation of IFN-I" = "#6FA7A0",
  "Inborn errors of cell death" = "#E9822A",
  "Immune metabolic" = "#E9C348", "Inflammasome IL-1β" = "#8A9BC8",
  "Arachidonic acid metabolism" = "#5AAFC7",
  "Endolysosomal nucleic acid sensing" = "#C05C86",
  "Ca2+ flux-PLC" = "#9DC6E1", "Osteoclast function" = "#F0B6C7",
  "NF-κB pathway" = "#E6AE70"
)

map_disease_name <- function(d) {
  d <- gsub("_", " ", d)
  d <- gsub("NFKB pathway", "NF-κB pathway", d)
  d <- gsub("TBK1 IRF3", "TBK1_IRF3", d)
  return(d)
}

col_lum <- function(hex) {
  rgb <- col2rgb(hex) / 255
  0.299 * rgb[1] + 0.587 * rgb[2] + 0.114 * rgb[3.5]
}

am_tab <- bidding_act$result %>%
  select(disease, module, activity) %>%
  mutate(mod_num = as.numeric(sub("M", "", module)),
         disease_color = disease_colors[map_disease_name(disease)],
         disease_color = ifelse(is.na(disease_color), "#FFFFFF", disease_color)) %>%
  arrange(mod_num, disease) %>%
  mutate(disease = factor(disease, levels = rev(disease)))

im_tab <- bidding_inh$result %>%
  select(disease, module, activity) %>%
  mutate(mod_num = as.numeric(sub("M", "", module)),
         disease_color = disease_colors[map_disease_name(disease)],
         disease_color = ifelse(is.na(disease_color), "#FFFFFF", disease_color)) %>%
  arrange(mod_num, disease) %>%
  mutate(disease = factor(disease, levels = rev(disease)))

plot_table <- function(df, title_text, header_fill = "#333333",
                       col1_w = 1, col2_w = 3, gap = 0.15) {
  n_rows <- nrow(df)
  y_vals <- as.numeric(df$disease)
  
  # column centers
  x1 <- col1_w / 2
  x2 <- col1_w + gap + col2_w / 2
  x_max <- col1_w + gap + col2_w
  
  col1_df <- data.frame(
    x = x1, y = y_vals, label = df$module,
    fill = "#FFFFFF", text_color = "black",
    tile_w = col1_w * 0.98, stringsAsFactors = FALSE
  )
  col2_df <- data.frame(
    x = x2, y = y_vals, label = as.character(df$disease),
    fill = df$disease_color,
    text_color = ifelse(sapply(df$disease_color, col_lum) > 0.5, "black", "white"),
    tile_w = col2_w * 0.98, stringsAsFactors = FALSE
  )
  plot_df <- rbind(col1_df, col2_df)
  
  ggplot(plot_df, aes(x = x, y = y)) +
    geom_tile(aes(fill = fill, width = tile_w, height = 0.98),
              color = "black", linewidth = 0.4) +
    scale_fill_identity() +
    geom_text(aes(label = label, color = text_color),
              size = 3, hjust = 0.5, vjust = 0.5) +
    scale_color_identity() +
    # header
    annotate("rect", xmin = 0.01, xmax = col1_w * 0.99,
             ymin = n_rows + 0.5, ymax = n_rows + 1.5,
             fill = header_fill, color = "black", linewidth = 0.4) +
    annotate("rect", xmin = col1_w + gap + 0.01, xmax = x_max - 0.01,
             ymin = n_rows + 0.5, ymax = n_rows + 1.5,
             fill = header_fill, color = "black", linewidth = 0.4) +
    annotate("text", x = x1, y = n_rows + 1, label = title_text,
             color = "white", size = 3.5, fontface = "bold") +
    annotate("text", x = x2, y = n_rows + 1, label = "Disease",
             color = "white", size = 3.5, fontface = "bold") +
    scale_x_continuous(limits = c(0, x_max), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, n_rows + 1.5), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(10, 10, 10, 10))
}

p_am <- plot_table(am_tab, "AM")
p_im <- plot_table(im_tab, "IM")

p_am + p_im + plot_layout(ncol = 2, widths = c(1, 1))


####writing result
##(1) Factor ~ cell type
labe218<-read.csv("E:\\AID cohort\\code\\NMF\\Spectra\\lam0001_dict242_218factors_labeled.csv")
library(tidyr)
labe218 <-labe218%>%
  extract(X, into = c("Factor", "Pathway"),
          regex = "^(Factor_\\d+)[_\\s]+(.+)$",
          remove = FALSE)
labe218<-labe218[,c("Factor","Pathway")]
library(openxlsx)
write.xlsx(labe218,file="E:\\AID cohort\\SupplememtTablexlsx\\Table 5_Annotated 218 factors.xlsx")
##(2) NMF result
# iterate over each module, collect unique factors from factor1+factor2, split into Factor and label
library(tidyr)
library(dplyr)

Act_result <- do.call(rbind, lapply(names(network_act_pw), function(mod) {
  edges   <- network_act_pw[[mod]]$edges_in_graph
  factors <- unique(c(as.character(edges$factor1), as.character(edges$factor2)))
  tibble(X = factors, Module = mod) %>%
    extract(X, into = c("Factor", "label"), regex = "^(Factor_\\d+)[_[:space:]]+(.+)$")
}))

Inh_result <- do.call(rbind, lapply(names(network_inh_pw), function(mod) {
  edges   <- network_inh_pw[[mod]]$edges_in_graph
  factors <- unique(c(as.character(edges$factor1), as.character(edges$factor2)))
  # strip leading Inhibit_ first
  factors_clean <- sub("^Inhibit_", "", factors)
  tibble(X = factors_clean, Module = mod) %>%
    extract(X, into = c("Factor", "label"), regex = "^(Factor_\\d+)[_[:space:]]+(.+)$")
}))

# join with bidding result using original M1‑M12, then add prefix A / I
Act_result <- Act_result %>%
  left_join(bidding_act$result, by = c("Module" = "module")) %>%
  mutate(Module = paste0("A", Module))
Inh_result <- Inh_result %>%
  left_join(bidding_inh$result, by = c("Module" = "module")) %>%
  mutate(Module = paste0("I", Module))

# check split failures
Act_bad <- Act_result %>% filter(is.na(Factor) | is.na(label))
Inh_bad <- Inh_result %>% filter(is.na(Factor) | is.na(label))

All_module_result <- bind_rows(Act_result, Inh_result)
write.xlsx(All_module_result,file="E:\\AID cohort\\SupplememtTablexlsx\\Table 8_NMF module activity.xlsx")
