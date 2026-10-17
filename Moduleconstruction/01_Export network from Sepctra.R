##### read spectra factors and gene dict to export  (pathway-gene-weight network,only contain genes in the given pathways)
library(dplyr)
library(readr)
library(lme4)
library(lmerTest)
library(stringr)
library(tidyr)
library(data.table)
library(openxlsx)
library(tibble)

setwd("E:\\AID cohort\\code\\NMF\\Spectra")

#--------1 export  network
#spectra_factors<-read.csv("lam0001_dict242_spectra_run_factors.csv")#258 factor * 21002 genes
#length(unique(spectra_factors$X))#258
mydict242<-read.xlsx("E:\\AID cohort\\code\\NMF\\Spectra\\mydict备份\\mydict_242.xlsx")
length(unique(mydict242$process))#242
spectra_factors_select<-read.csv("lam0001_dict242_218factors_labeled.csv")
length(unique(spectra_factors_select$X))#2218
#gsub("^Factor_\\d+_\\s*", "", "Factor_67_pDC_Type_I_interferon_signaling" )

spectra_factors_select$Y<-sapply(spectra_factors_select$X,function(x){
  trimws(gsub("^Factor_\\d+[_\\s]*", "", x))
})
sum(is.na(match(spectra_factors_select$Y,mydict242$process)))# 0 completing match

Spec_network<-data.frame()
for(i in 1:nrow(spectra_factors_select)){
  paths<-spectra_factors_select$Y[i]
  genes<-filter(mydict242,process==paths)%>%dplyr::select(gene)%>%pull()
  resu<-spectra_factors_select[i,na.omit(match(genes,colnames(spectra_factors_select)))]
  resu2<-data.frame(source=rep(spectra_factors_select$X[i],ncol(resu)),
                    target=names(resu),
                    weight=as.numeric(resu))
  Spec_network<-rbind(Spec_network,resu2)
}
dim(Spec_network)#15139
length(unique(Spec_network$source))#218
saveRDS(Spec_network,file="E:\\AID cohort\\code\\RegVelo\\Spect_network.RDS")


#--------2 produce wide format LMEM for NMF 1017 revised
library(data.table)
#file<-fread("E:\\AID cohort\\code\\NMF\\0908New_condition\\lam0001_dict242_lmem_results_annotated_6groups.txt")
file<-fread("E:\\AID cohort\\code\\NMF\\Spectra\\LMEM建模结果备份\\1016lam0001_dict242_lmem_results_annotated.txt")
Lmem_long<-as.data.frame(file)
head(Lmem_long)
sum(Lmem_long$fdr_global>0.05)#1379, 20261016(2470)
Lmem_long<-filter(Lmem_long,fdr_global<=0.05,celltype!="all")
class(Lmem_long$fdr_global)
min(Lmem_long$fdr_global)
length(unique(Lmem_long$factor))#144 118 20261016(142)
unique(Lmem_long$disease)#16 9 20261016(16)

Lmem_wide <- Lmem_long %>%
  dplyr::select(disease, factor, beta) %>%
  pivot_wider(names_from = disease, values_from = beta, values_fill = 0) %>%
  column_to_rownames("factor")
dim(Lmem_wide)
sum(Lmem_wide<0)# 299,  20261016(280) existing negative beta value  could turn to doublee matrix, proposed by glioma(cancer cell)

new_Lmem_wide<-data.frame()
for(i in 1:nrow(Lmem_wide)){
  values<-Lmem_wide[i,];cut_offmin<-min(values)
  if(cut_offmin>=0){
    print(paste0("-----",rownames(Lmem_wide)[i],"  all positive beta-----"))
  res<-values
  }else{
    value1<-as.numeric(as.matrix(values))
    values2<-ifelse(value1<0,0,value1)
    values3 <- ifelse(value1<0, abs(value1),0)#inhibit
    res<-as.data.frame(rbind(values2,values3))
    colnames(res)<-colnames(values)
    rownames(res)<-c(
      rownames(values),
      paste0("Inhibit_",rownames(values))
    )
    print(paste0("-----",rownames(Lmem_wide)[i],"turn to negative and positive both existing"))
  }
  new_Lmem_wide<-rbind(new_Lmem_wide,res)
}
dim(new_Lmem_wide)#230  #175  20261016(225)


#save(new_Lmem_wide,file="E:\\AID cohort\\code\\NMF\\0908New_condition\\sigLeme_widebeta.RDS")
save(new_Lmem_wide,file="E:\\AID cohort\\code\\NMF\\sigLeme_widebeta.RDS")



