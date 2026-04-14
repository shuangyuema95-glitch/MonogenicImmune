library(tidyverse)

multiplesheets <- function(fname) {
  
  # getting info about all excel sheets
  sheets <- readxl::excel_sheets(fname)
  tibble <- lapply(sheets, function(x) readxl::read_excel(fname, sheet = x))
  data_frame <- lapply(tibble, as.data.frame)
  
  # assigning names to data frames
  names(data_frame) <- sheets
  
  # print data frame
  print(data_frame)
}

Clist<-multiplesheets("E:\\AID cohort\\data\\sigenature collection.xlsx")

#(1) Nano Straing  human immunological panel V2
x1<-Clist[[1]]
Nano<-data.frame()
for(i in 2:ncol(x1)){
  x1[,c(1,i)]->data
  data2<-data.frame(gene=data[,1][which(data[,2]=="+")])
  data2$pathway<-rep(colnames(x1)[i],nrow(data2))
  Nano<-rbind(Nano,data2)
}
Nano$source<-rep("NanoString Human immunology panel 2")
dim(Nano)
unique(Nano$pathway)

#(2) PMID: 27568558 
Ur2016<-Clist[[2]]

#(3) Nature Medicine, PMID: 41526507
NM<-Clist[[3]]
NMp<-data.frame()
for(i in 2:ncol(NM)){
  NM[,c(1,i)]->data
  data2<-data.frame(gene=data[,1][which(data[,2]=="TRUE")])
  data2$pathway<-rep(colnames(NM)[i],nrow(data2))
  NMp<-rbind(NMp,data2)
}
NMp$source<-rep("Jiménez-Gracia et al.")

#(4) Nature Cancer, PMID: 40854985
NCC<-Clist[[4]]
NCCp<-data.frame()
for(i in 1:ncol(NCCp)){
  NCCp<-rbind(NCCp,data.frame(gene=as.character(na.omit(NCC[,i])),
             pathway=colnames(NCC)[i]))}
NCCp$source<-rep("Han et al.")

Y<-do.call("rbind",list(Nano,Ur2016,NMp,NCCp))
unique(Y$pathway)


#(5) got gene list form MsigDB datdabase
library(msigdbr)
species <- "Homo sapiens"
msig <- msigdbr(species = species)

# Hallmark (H)
hallmark_sets <- split(
  msig$gene_symbol[msig$gs_cat == "H"], 
  msig$gs_name[msig$gs_cat == "H"]
)


# KEGG (C2:CP:KEGG)
kegg_sets <- split(
  msig$gene_symbol[msig$gs_subcat == "CP:KEGG"], 
  msig$gs_name[msig$gs_subcat == "CP:KEGG"]
)

# Reactome (C2:CP:REACTOME)
reactome_sets <- split(
  msig$gene_symbol[msig$gs_subcat == "CP:REACTOME"], 
  msig$gs_name[msig$gs_subcat == "CP:REACTOME"]
)

# GO BP (C5:GO:BP)
go_sets <- split(
  msig$gene_symbol[msig$gs_subcat == "GO:BP"], 
  msig$gs_name[msig$gs_subcat == "GO:BP"]
)

# Immunologic signatures (C7)
c7_sets <- split(
  msig$gene_symbol[msig$gs_cat == "C7"], 
  msig$gs_name[msig$gs_cat == "C7"]
)

hallmark_sets <- lapply(hallmark_sets, unique)
kegg_sets     <- lapply(kegg_sets, unique)
reactome_sets <- lapply(reactome_sets, unique)
go_sets       <- lapply(go_sets, unique)
c7_sets       <- lapply(c7_sets, unique)


all_sets <- c(
  hallmark_sets,
  kegg_sets,
  reactome_sets,
  go_sets,
  c7_sets
)



convert_to_dataframe <- function(sets_list) {
  #  purrr::imap ,x is the vector value of each element , y is the name of each element
  all_genes_df <- purrr::imap_dfr(sets_list, ~tibble(
    gene = .x, 
    pathway = .y,
    source = strsplit(.y, "_")[[1]][1]  
  ))
  
  return(all_genes_df)
}

all_genes_df <- convert_to_dataframe(all_sets)
head(all_genes_df)


save(all_genes_df,file="E:\\AID cohort\\code\\result\\genelistDF.Rdata")


