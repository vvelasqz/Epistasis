#define main operations to calculate epistasis

time=c(0, 16, 20, 24, 32, 48)
gen=c("wt", "mla6", "rar3","bln1", "dm")
combi <-combn(c("wt", "mla6", "rar3","bln1", "dm"),2)
comb<- apply(combi, 2, FUN=function(x)paste0(x[1],"_", x[2]))
FC_padj_table<- read.csv("~/iowa_state/lab/RNAseq/DESeq2 HvR3/hv_R3_genes_deseq2_results_pairwise_genotype_logfc_padj_tax_sp.csv")
#FC_padj_table<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/hv_R2_genes_deseq2_results_pairwise_genotype_logfc_padj_tax_sp.csv")
#FC_padj_table<-read.csv("~/iowa_state/lab/RNAseq/DESeq2 analysis HvV2/hv_V2_genes_deseq2_results_pairwise_genotype_FC_padj_tax_sp.csv", stringsAsFactors = F)
FC_padj_table[,grep("FoldChange", colnames(FC_padj_table))][is.na(FC_padj_table[,grep("FoldChange", colnames(FC_padj_table))])]<- 0
FC_padj_table[,grep("padj", colnames(FC_padj_table))][is.na(FC_padj_table[,grep("padj", colnames(FC_padj_table))])]<- 1

#create a function that compares using the 2 genotypes, and p value, report with FC and fix by num and den, and if DE or !DE

DE_list<- function(DE_table, num, den, padj, DE=T, FC=0, comb, times){
  library(textclean)
  if(paste0(num, "_", den) %in% comb){
    DE_table_mod<- DE_table[, c(1, grep(paste(paste0(num, "_", den), collapse = "|"), colnames(DE_table)))]
  }else{
    DE_table_mod<- DE_table[, c(1, grep(paste(paste0(den, "_", num), collapse = "|"), colnames(DE_table)))]
    DE_table_mod[, grep("FoldChange", colnames(DE_table_mod))]<- -DE_table_mod[, grep("FoldChange", colnames(DE_table_mod))]
    colnames(DE_table_mod)<- swap(colnames(DE_table_mod), num, den)
  }
  #DE_table_mod<- DE_table_mod[complete.cases(DE_table_mod),]
  result<- list()
  for(t in times){
    if(DE==T){
      if(FC>0){
        result[[paste0("t",t)]]<- DE_table_mod[DE_table_mod[,paste0("t",t,"_", num, "_", den,"_padj")]< padj &
                                                 DE_table_mod[,paste0("t",t,"_", num, "_", den,"_log2FoldChange")] > 0, c(1, grep(paste0("t",t,"_", num, "_", den), colnames(DE_table_mod)))]
      }else if(FC<0){
        result[[paste0("t",t)]]<- DE_table_mod[DE_table_mod[,paste0("t",t,"_", num, "_", den,"_padj")]< padj &
                                                 DE_table_mod[,paste0("t",t,"_", num, "_", den,"_log2FoldChange")] < 0, c(1, grep(paste0("t",t,"_", num, "_", den), colnames(DE_table_mod)))]
      }else{
        result[[paste0("t",t)]]<- DE_table_mod[DE_table_mod[,paste0("t",t,"_", num, "_", den,"_padj")]< padj, c(1, grep(paste0("t",t,"_", num, "_", den), colnames(DE_table_mod)))]
      }
    }else{
      result[[paste0("t",t)]]<- DE_table_mod[DE_table_mod[,paste0("t",t,"_", num, "_", den,"_padj")]> padj, c(1, grep(paste0("t",t,"_", num, "_", den), colnames(DE_table_mod)))]
      #result[[paste0("t",t)]]<- result[[paste0("t",t)]][!duplicated(result[[paste0("t",t)]]$gene),]
    }
    #result[[paste0("t",t)]]<- result[[paste0("t",t)]][!is.na(result[[paste0("t",t)]]$gene),]
  }
  return(result)
}

#z<- result[["t20"]]

set_operations<- function(list1, list2, operation="intersection"){
  result<- list()
  if(operation=="intersection"){
    for(t in names(list1)){
      result[[t]]<- list1[[t]][list1[[t]]$gene %in% list2[[t]]$gene, ]
      if(sum(!colnames(list2[[t]]) %in% colnames(list1[[t]]))>0){
        result[[t]]<- cbind(result[[t]], list2[[t]][match(result[[t]]$gene, list2[[t]]$gene),-grep(paste(colnames(list1[[t]]), collapse = "|"), colnames(list2[[t]]))])
      }
    }
  }else if(operation=="minus"){
    for(t in names(list1)){
      result[[t]]<- list1[[t]][!list1[[t]]$gene %in% list2[[t]]$gene, ]
    }
  }else if(operation=="union"){
    library(dplyr)
    for(t in names(list1)){
      common_genes<- list1[[t]][list1[[t]]$gene %in% list2[[t]]$gene, "gene"]
      result[[t]]<- bind_rows(merge(list1[[t]][list1[[t]]$gene %in% common_genes, ], list2[[t]][match(common_genes, list2[[t]]$gene),], by=colnames(list1[[t]])[colnames(list1[[t]]) %in% colnames(list2[[t]])]),
                          list1[[t]][!list1[[t]]$gene %in% common_genes, ], list2[[t]][!list2[[t]]$gene %in% common_genes, ])
    }
  }
  return(result)
}

get_epistatic<-function(triple_union_list, thr=1, thr_prop=0.2, mode="dev"){
  result<- triple_union_list
  result_notEpi<-list()
  for(t in names(result)){
    result[[t]][,grep("FoldChange", colnames(result[[t]]))][is.na(result[[t]][,grep("FoldChange", colnames(result[[t]]))])]<- 0
    result[[t]]$expected<-result[[t]][,4]+result[[t]][,6]
    result[[t]]$deviation<-result[[t]][,2]-result[[t]][,4]-result[[t]][,6]
    result[[t]]$dev_prop<-abs(result[[t]]$deviation/(result[[t]][,2]+0.1))
    if(mode=="dev"){
      result_notEpi[[t]]<- result[[t]][abs(result[[t]]$deviation)<thr,]
      result[[t]]<- result[[t]][abs(result[[t]]$deviation)>thr,]
    }else{
      result_notEpi[[t]]<- result[[t]][result[[t]]$dev_prop<thr_prop,]
      result[[t]]<- result[[t]][result[[t]]$dev_prop>thr_prop,]
    }
  }
  return(list(result,result_notEpi))
}

#for barley
horvu_DE_table<- FC_padj_table[grep("HORVU",FC_padj_table$gene), ]
#Mla6 and Rar3
wt_mla6_DE<- DE_list(DE_table=horvu_DE_table, num="wt", den="mla6", padj=0.001, DE=T, comb=comb, times=time)
wt_rar3_DE<- DE_list(DE_table=horvu_DE_table, num="wt", den="rar3", padj=0.001, DE=T, comb=comb, times=time)
mla6_rar3_notDE<- DE_list(DE_table=horvu_DE_table, num="mla6", den="rar3", padj=0.001, DE=F, comb=comb, times=time)
mla6_rar3_DE_FC_plus<-DE_list(DE_table=horvu_DE_table, num="mla6", den="rar3", padj=0.001, DE=T, FC=1, comb=comb, times=time)
mla6_rar3_DE_FC_minus<-DE_list(DE_table=horvu_DE_table, num="mla6", den="rar3", padj=0.001, DE=T, FC=-1, comb=comb, times=time)

shared_eff_Mla6_Rar3<- set_operations(list1=wt_mla6_DE, list2=wt_rar3_DE, operation="intersection")
unique_eff_Mla6<- set_operations(list1=wt_mla6_DE, list2=shared_eff_Mla6_Rar3, operation="minus")
unique_eff_Rar3<- set_operations(list1=wt_rar3_DE, list2=shared_eff_Mla6_Rar3, operation="minus")

equal_eff_Mla6_Rar3<- set_operations(list1=shared_eff_Mla6_Rar3, list2=mla6_rar3_notDE, operation="intersection")
predominant_Mla6<- set_operations(list1=shared_eff_Mla6_Rar3, list2=mla6_rar3_DE_FC_plus, operation="intersection")
predominant_Rar3<- set_operations(list1=shared_eff_Mla6_Rar3, list2=mla6_rar3_DE_FC_minus, operation="intersection")


#Mla6 and Bln1
wt_mla6_DE<- DE_list(DE_table=horvu_DE_table, num="wt", den="mla6", padj=0.001, DE=T, comb=comb, times=time)
wt_bln1_DE<- DE_list(DE_table=horvu_DE_table, num="wt", den="bln1", padj=0.001, DE=T, comb=comb, times=time)
wt_dm_DE<- DE_list(DE_table=horvu_DE_table, num="wt", den="dm", padj=0.001, DE=T, comb=comb, times=time)

mla6_bln1_notDE<- DE_list(DE_table=horvu_DE_table, num="mla6", den="bln1", padj=0.001, DE=F, comb=comb, times=time)
mla6_dm_notDE<- DE_list(DE_table=horvu_DE_table, num="mla6", den="dm", padj=0.001, DE=F, comb=comb, times=time)
bln1_dm_notDE<- DE_list(DE_table=horvu_DE_table, num="bln1", den="dm", padj=0.001, DE=F, comb=comb, times=time)

#shared_eff_Mla6_Bln1<- set_operations(list1=wt_mla6_DE, list2=wt_bln1_DE, operation="intersection")
#shared_eff_Mla6_Bln1<- set_operations(list1=shared_eff_Mla6_Bln1, list2=mla6_bln1_notDE, operation="intersection")
epi_Mla6_Bln1_DE<-set_operations(list1=wt_dm_DE, list2=wt_mla6_DE, operation="union")
epi_Mla6_Bln1_DE<-set_operations(list1=epi_Mla6_Bln1_DE, list2=wt_bln1_DE, operation="union")
#change this so only er look at DE genes in wt dm comparisson
#epi_Mla6_Bln1<-set_operations(list1=epi_Mla6_Bln1, list2=wt_bln1_DE, operation="intersection")

epi_Mla6_Bln1_all_prop<- get_epistatic(epi_Mla6_Bln1_DE, thr_prop = 0.5, mode = "prop")
a1<-epi_Mla6_Bln1_all_prop[[1]][["t16"]]
epi_Mla6_Bln1_prop<-epi_Mla6_Bln1_all_prop[[1]]
additive_Mla6_Bln1_prop<-epi_Mla6_Bln1_all_prop[[2]]

epi_Mla6_Bln1_all_dev<- get_epistatic(epi_Mla6_Bln1_DE, thr = 1)
epi_Mla6_Bln1<-epi_Mla6_Bln1_all_dev[[1]]
additive_Mla6_Bln1<-epi_Mla6_Bln1_all_dev[[2]]


epi_Mla6_Bln1<-epi_Mla6_Bln1_all_prop[[1]]
additive_Mla6_Bln1<-epi_Mla6_Bln1_all_prop[[2]]

a<-epi_Mla6_Bln1[["t16"]]
b<-additive_Mla6_Bln1[["t20"]]

mla6_dm_DE_FC_plus<- DE_list(DE_table=horvu_DE_table, num="mla6", den="dm", padj=0.001, DE=T, FC=1, comb=comb, times=time)
mla6_dm_DE_FC_minus<- DE_list(DE_table=horvu_DE_table, num="mla6", den="dm", padj=0.001, DE=T, FC=-1, comb=comb, times=time)

a1<-mla6_dm_DE_FC_plus[["t20"]]
b1<-mla6_dm_DE_FC_minus[["t20"]]

bln1_dm_DE_FC_plus<- DE_list(DE_table=horvu_DE_table, num="bln1", den="dm", padj=0.001, DE=T, FC=1, comb=comb, times=time)
bln1_dm_DE_FC_minus<- DE_list(DE_table=horvu_DE_table, num="bln1", den="dm", padj=0.001, DE=T, FC=-1, comb=comb, times=time)

positive<- set_operations(list1=epi_Mla6_Bln1, list2=bln1_dm_DE_FC_minus, operation="intersection")
positive<- set_operations(list1=positive, list2=mla6_dm_DE_FC_minus, operation="intersection")

apos<-positive[["t16"]]

pseudo_masked<- set_operations(list1=epi_Mla6_Bln1, list2=bln1_dm_DE_FC_plus, operation="intersection")
pseudo_masked<- set_operations(list1=pseudo_masked, list2=mla6_dm_DE_FC_minus, operation="intersection")

pseudo_masked2<- set_operations(list1=epi_Mla6_Bln1, list2=bln1_dm_DE_FC_minus, operation="intersection")
pseudo_masked2<- set_operations(list1=pseudo_masked2, list2=mla6_dm_DE_FC_plus, operation="intersection")
pseudo_masked<- set_operations(list1=pseudo_masked, list2=pseudo_masked2, operation="union")


d<-bln1_dm_DE_FC_plus[["t20"]]
bPM<-pseudo_masked[["t20"]]

negative<- set_operations(list1=epi_Mla6_Bln1, list2=bln1_dm_DE_FC_plus, operation="intersection")
negative<- set_operations(list1=negative, list2=mla6_dm_DE_FC_plus, operation="intersection")

c<-negative[["t20"]]

#masked<- set_operations(list1=wt_dm_DE, list2=wt_mla6_DE, operation="intersection")
#masked<- set_operations(list1=masked, list2=shared_eff_Mla6_Bln1, operation="minus")
masked<- set_operations(list1=epi_Mla6_Bln1, list2=mla6_dm_notDE, operation="intersection")

e<-masked[["t20"]]

#suppression<- set_operations(list1=wt_bln1_DE, list2=wt_dm_DE, operation="intersection")
#suppression<- set_operations(list1=suppression, list2=shared_eff_Mla6_Bln1, operation="minus")
suppression<- set_operations(list1=epi_Mla6_Bln1, list2=bln1_dm_notDE, operation="intersection")

f<-suppression[["t20"]]

sum(e$gene %in% f$gene)
g<-e$gene[e$gene %in% f$gene]

#If we want to separate further symmetric interaction and remove those from masked and suppression
symmetric<- set_operations(list1=masked, list2=suppression, operation="intersection")
masked<- set_operations(list1=masked, list2=symmetric, operation="minus")
suppression<- set_operations(list1=suppression, list2=symmetric, operation="minus")

#intersections
#make a function to classify the genes with the epistasis patters by timepoint
#find if one gene has more than one pattern

combine_lists<- function(oldlist){
  library(data.table)
  new_list<- list()
  for (name in names(oldlist)) {
    if(nrow(oldlist[[name]]>0)){
      new_list[[name]]<- data.frame(gene=oldlist[[name]]$gene, time=name, epistasis=deparse(substitute(oldlist)))
    }
  }
  combined<- rbindlist(new_list)
  return(combined)
}

#shared_Mla6_Rar3<- combine_lists(shared_eff_Mla6_Rar3)
unique_Mla6<- combine_lists(unique_eff_Mla6)
unique_Rar3<- combine_lists(unique_eff_Rar3)
equal_Mla6_Rar3<- combine_lists(equal_eff_Mla6_Rar3)
predom_Mla6<- combine_lists(predominant_Mla6)
predom_Rar3<- combine_lists(predominant_Rar3)

library(tidyverse)
#combined_tables_Mla6_Rar3<- rbind(shared_Mla6_Rar3,unique_Mla6, unique_Rar3, equal_Mla6_Rar3, predom_Mla6, predom_Rar3)
combined_tables_Mla6_Rar3<- rbind(unique_Mla6, unique_Rar3, equal_Mla6_Rar3, predom_Mla6, predom_Rar3)
#combined_tables_Mla6_Rar3<- combined_tables_Mla6_Rar3 %>% group_by(gene, time) %>% summarise(epistasis = paste0(epistasis, collapse = ";"))
spread_tables_Mla6_Rar3<- spread(combined_tables_Mla6_Rar3, key = time, value = epistasis)
summary_barley_Mla6_Rar3<- t(rbindlist(sapply(X = spread_tables_Mla6_Rar3[,-1], FUN=function(x)data.frame(as.list(table(x)))), fill = TRUE))
colnames(summary_barley_Mla6_Rar3)<- paste0("t", time)
#summary_barley_Mla6_Rar3<- sapply(X = spread_tables_Mla6_Rar3[,-1], FUN = table)


addi_Mla6_Bln1<- combine_lists(additive_Mla6_Bln1)
pos<- combine_lists(positive)
pseudo<- combine_lists(pseudo_masked)
neg<- combine_lists(negative)
masked_epi<- combine_lists(masked)
supp<- combine_lists(suppression)
symmetric_Mla6_Bln1<- combine_lists(symmetric)

#combined_tables_Mla6_Bln1<- rbind(addi_Mla6_Bln1, pos, pseudo, neg, masked_epi, supp)
combined_tables_Mla6_Bln1<- rbind(addi_Mla6_Bln1, pos, pseudo, neg, masked_epi, supp, symmetric_Mla6_Bln1)
combined_tables_Mla6_Bln1<- combined_tables_Mla6_Bln1 %>% group_by(gene, time) %>% summarise(epistasis = paste0(epistasis, collapse = ";"))
spread_tables_Mla6_Bln1<- spread(combined_tables_Mla6_Bln1, key = time, value = epistasis)

summary_barley_Mla6_Bln1<- sapply(X = spread_tables_Mla6_Bln1[,-1], FUN = table)


#spread_tables_Mla6_Bln1[spread_tables_Mla6_Bln1=="additive_Mla6_Bln1"]<-""
spread_tables_Mla6_Bln1$consensus<- apply(spread_tables_Mla6_Bln1[,3:7], 1, 
                                          FUN=function(x){y=table(unlist(x))
                                          consensus=names(y)[which(y==max(y))]
                                          ifelse("additive_Mla6_Bln1" %in% consensus & length(consensus)>1,
                                                 paste(consensus[!consensus=="additive_Mla6_Bln1"], collapse = ";"),
                                                 paste(consensus, collapse = ";"))})

unique(spread_tables_Mla6_Bln1$consensus)
spread_tables_Mla6_Bln1<- spread_tables_Mla6_Bln1[spread_tables_Mla6_Bln1$consensus!="",]
table(spread_tables_Mla6_Bln1$consensus)
#spread_tables_Mla6_Bln1[spread_tables_Mla6_Bln1$consensus=="", "consensus"]<- "additive_Mla6_Bln1"
#spread_tables_Mla6_Bln1[is.na(spread_tables_Mla6_Bln1)]<- "additive_Mla6_Bln1"


spread_tables_Mla6_Rar3$consensus<- apply(spread_tables_Mla6_Rar3[,3:7], 1, 
                                          FUN=function(x){y=table(unlist(x)) 
                                          paste(names(y)[which(y==max(y))], collapse = ";")})


#save the lists
#deparse(substitute(inter_common_epi_mla6_Bln1))
save_lists<- function(list, path){
  for (name in names(list)) {
    if(nrow(list[[name]]>0)){
    write.csv(list[[name]], paste0(path, name,"_",deparse(substitute(list)),".csv"), row.names = F)
    }}
  saveRDS(list, paste0(path, deparse(substitute(list)),".RDS"))
}

#save_lists(list=shared_eff_Mla6_Rar3, path="~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/")
# save_lists(list=unique_eff_Mla6, path="~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/")
# save_lists(list=unique_eff_Rar3, path="~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/")
# save_lists(list=equal_eff_Mla6_Rar3, path="~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/")
# save_lists(list=predominant_Mla6, path="~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/")
# save_lists(list=predominant_Rar3, path="~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/")
# write.csv(spread_tables_Mla6_Rar3, "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/spread_tables_Mla6_Rar3.csv")
# 
# save_lists(list=epi_Mla6_Bln1, path="~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/")
# save_lists(list=additive_Mla6_Bln1, path="~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/")
# 
# save_lists(list=positive, path="~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/")
# save_lists(list=pseudo_masked, path="~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/")
# save_lists(list=negative, path="~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/")
# save_lists(list=masked, path="~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/")
# save_lists(list=suppression, path="~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/")
# save_lists(list=symmetric, path="~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/")
# write.csv(spread_tables_Mla6_Bln1, "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/spread_tables_Mla6_Bln1.csv")
# write.csv(summary_barley_Mla6_Bln1, "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/summary_barley_Mla6_Bln1.csv")
# write.csv(summary_barley_Mla6_Rar3, "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/summary_barley_Mla6_Rar3.csv")

# save_lists(list=mla6_dm_DE_FC_plus, path="~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/")
# save_lists(list=mla6_dm_DE_FC_minus, path="~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/")
# save_lists(list=bln1_dm_DE_FC_plus, path="~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/")
# save_lists(list=bln1_dm_DE_FC_minus, path="~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/")

#HVR3
save_lists(list=unique_eff_Mla6, path="~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/")
save_lists(list=unique_eff_Rar3, path="~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/")
save_lists(list=equal_eff_Mla6_Rar3, path="~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/")
save_lists(list=predominant_Mla6, path="~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/")
save_lists(list=predominant_Rar3, path="~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/")
write.csv(spread_tables_Mla6_Rar3, "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/spread_tables_Mla6_Rar3.csv")
write.csv(summary_barley_Mla6_Rar3, "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/summary_barley_Mla6_Rar3.csv")

save_lists(list=epi_Mla6_Bln1, path="~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/")
save_lists(list=additive_Mla6_Bln1, path="~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/")

save_lists(list=positive, path="~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/")
save_lists(list=pseudo_masked, path="~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/")
save_lists(list=negative, path="~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/")
save_lists(list=masked, path="~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/")
save_lists(list=suppression, path="~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/")
save_lists(list=symmetric, path="~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/")
write.csv(spread_tables_Mla6_Bln1, "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/spread_tables_Mla6_Bln1.csv")
write.csv(summary_barley_Mla6_Bln1, "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/summary_barley_Mla6_Bln1.csv")
