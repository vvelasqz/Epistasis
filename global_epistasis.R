library(reshape2)
library(tidyverse)
# global epistasis calculations
#fc_table<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/hv_R2_genes_deseq2_results_pairwise_genotype_logfc_tax_sp.csv")
fc_table<- read.csv("~/iowa_state/lab/RNAseq/DESeq2 HvR3/hv_R3_genes_deseq2_results_pairwise_genotype_logfc_tax_sp.csv")
fc_table_mla6_bln1<- melt(fc_table[,grepl("gene|wt_dm", colnames(fc_table))])
fc_table_mla6<- melt(fc_table[,grepl("gene|wt_mla6", colnames(fc_table))])
fc_table_bln1<- melt(fc_table[,grepl("gene|wt_bln1", colnames(fc_table))])
#melt(fc_table_mla6_bln1)
time=c(0, 16, 20, 24, 32, 48)
gen=c("wt", "mla6", "rar3","bln1", "dm")

pred_obs<- data.frame(gene=rep(fc_table$gene, length(time)), 
                      time=rep(paste0("t",time), each=nrow(fc_table)))
pred_obs$org<- ifelse(grepl("BLGH", pred_obs$gene, ignore.case=TRUE), "Bgh", "Barley")
pred_obs$pred_addi<- fc_table_mla6$value + fc_table_bln1$value
pred_obs$obs<-fc_table_mla6_bln1$value
pred_obs$deviation<- pred_obs$obs - pred_obs$pred_addi
pred_obs<- pred_obs[!is.na(pred_obs$deviation),]

correlations<- pred_obs %>% group_by(org, time) %>% summarise(cor=cor(pred_addi, deviation))
correlations$cor<- round(correlations$cor, digits = 4)

correlations_all<- pred_obs %>% group_by(org) %>% summarise(cor=cor(pred_addi, deviation))
correlations_all$time<-"all"
correlations_all$cor<- round(correlations_all$cor, digits = 4)
correlations_all<-correlations_all[,c(1,3,2)]
correlations<- rbind(correlations, correlations_all)


#filter only DE genes and calculate slope of the line bgh 0.005, barley 0.001
#deg_table<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/hv_R2_genes_deseq2_results_pairwise_genotype_padj_tax_sp.csv")
deg_table<- read.csv("~/iowa_state/lab/RNAseq/DESeq2 HvR3/hv_R3_genes_deseq2_results_pairwise_genotype_padj_tax_sp.csv")
deg_table<- deg_table[,grepl("gene|wt_dm|wt_mla6|wt_bln1", colnames(deg_table))]

deg_table_mla6_bln1<- melt(deg_table[,grepl("gene|wt_dm", colnames(deg_table))])
deg_table_mla6<- melt(deg_table[,grepl("gene|wt_mla6", colnames(deg_table))])
deg_table_bln1<- melt(deg_table[,grepl("gene|wt_bln1", colnames(deg_table))])

deg<- data.frame(gene=rep(deg_table$gene, length(time)),
                 time=rep(paste0("t",time), each=nrow(deg_table)),
                 padj_mla6=deg_table_mla6$value, padj_bln1=deg_table_bln1$value,
                 padj_dm=deg_table_mla6_bln1$value)
deg$org<- ifelse(grepl("BLGH", deg$gene, ignore.case=TRUE), "Bgh", "Barley")
deg<- deg[complete.cases(deg),]
deg$de<- ifelse(deg$org=="Bgh", ifelse(deg$padj_mla6<0.005 | deg$padj_bln1<0.005 | deg$padj_dm<0.005, T, F),
                ifelse(deg$padj_mla6<0.001 | deg$padj_bln1<0.001 | deg$padj_dm<0.001, T, F))

pred_obs_deg<- merge(pred_obs, deg[, c(1,2,7)], by=c("gene", "time"))
  
pdf("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/global_epistasis/global_epistasis_deGenes.pdf", width = 18, height = 12, fonts = "ArialMT", pointsize = 18)

ggplot(data = pred_obs_deg[pred_obs_deg$de,], aes_string(x = "pred_addi", y = "deviation", color = "time")) + theme_bw(base_size=40) +
  geom_point(size = 2) + coord_fixed()+ scale_color_brewer(palette = "Dark2") +facet_wrap(~time)

ggplot(data = pred_obs_deg[pred_obs_deg$de,], aes_string(x = "pred_addi", y = "deviation", color = "time")) + theme_bw(base_size=40) +
  geom_point(size = 2) + coord_fixed()+ scale_color_brewer(palette = "Dark2") +facet_wrap(~org) +
  geom_smooth(method = "lm", se = FALSE)

ggplot(data = pred_obs_deg[pred_obs_deg$de & pred_obs_deg$org=="Barley",], aes_string(x = "pred_addi", y = "deviation", color = "time")) + theme_bw(base_size=40) +
  geom_point(size = 2) + coord_fixed()+ scale_color_brewer(palette = "Dark2") 

ggplot(data = pred_obs_deg[pred_obs_deg$de & pred_obs_deg$org=="Barley",], aes_string(x = "pred_addi", y = "deviation", color = "time")) + theme_bw(base_size=40) +
  geom_point(size = 2) + coord_fixed()+ scale_color_brewer(palette = "Dark2") +facet_wrap(~time)

ggplot(data = pred_obs_deg[pred_obs_deg$de & pred_obs_deg$org=="Bgh",], aes_string(x = "pred_addi", y = "deviation", color = "time")) + theme_bw(base_size=40) +
  geom_point(size = 2) + coord_fixed()+ scale_color_brewer(palette = "Dark2")

ggplot(data = pred_obs_deg[pred_obs_deg$de & pred_obs_deg$org=="Bgh",], aes_string(x = "pred_addi", y = "deviation", color = "time")) + theme_bw(base_size=40) +
  geom_point(size = 2) + coord_fixed()+ scale_color_brewer(palette = "Dark2") +facet_wrap(~time)

dev.off()

global_epistasis<- pred_obs_deg[pred_obs_deg$de,] %>% group_by(org, time) %>% 
  summarise(slope=lm(deviation ~ pred_addi)$coefficients[2])
global_epistasis$slope<- round(global_epistasis$slope, digits = 4)

#model<- lm(pred_obs_deg$deviation ~ pred_obs_deg$pred_addi)
#model$coefficients[2]

global_epistasis_all<- pred_obs_deg[pred_obs_deg$de,] %>% group_by(org) %>% 
  summarise(slope=lm(deviation ~ pred_addi)$coefficients[2])
global_epistasis_all$time<-"all"
global_epistasis_all$slope<- round(global_epistasis_all$slope, digits = 4)
global_epistasis_all<-global_epistasis_all[,c(1,3,2)]
global_epistasis<- rbind(global_epistasis, global_epistasis_all)

write.csv(global_epistasis, "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/global_epistasis/global_epsitasis_DEgenes.csv")



#filter only epistatic genes and calculate slope of the line 

#deg_table<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/hv_R2_genes_deseq2_results_pairwise_genotype_padj_tax_sp.csv")
spread_tables_Mla6_Bln1<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/spread_tables_Mla6_Bln1.csv", row.names = 1)
# spread_tables_Mla6_Bln1[spread_tables_Mla6_Bln1=="additive_Mla6_Bln1"]<-NA
# spread_tables_Mla6_Bln1$consensus<- apply(spread_tables_Mla6_Bln1[,2:7], 1, 
#                                           FUN=function(x){y=table(unlist(x)) 
#                                           paste(names(y)[which(y==max(y))], collapse = ";")})
# 
# spread_tables_Mla6_Bln1[spread_tables_Mla6_Bln1$consensus=="", "consensus"]<- "additive_Mla6_Bln1"
# spread_tables_Mla6_Bln1[is.na(spread_tables_Mla6_Bln1)]<- "additive_Mla6_Bln1"
spread_tables_Mla6_Bln1<- spread_tables_Mla6_Bln1[spread_tables_Mla6_Bln1$consensus!="additive_Mla6_Bln1",]

bgh_spread_tables_Mla6_Bln1<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/bgh/HVR3/bgh_spread_tables_Mla6_Bln1.csv", row.names = 1)
# bgh_spread_tables_Mla6_Bln1[bgh_spread_tables_Mla6_Bln1=="additive_Mla6_Bln1"]<-NA
# bgh_spread_tables_Mla6_Bln1$consensus<- apply(bgh_spread_tables_Mla6_Bln1[,2:7], 1, 
#                                               FUN=function(x){y=table(unlist(x)) 
#                                               paste(names(y)[which(y==max(y))], collapse = ";")})
# 
# bgh_spread_tables_Mla6_Bln1[bgh_spread_tables_Mla6_Bln1$consensus=="", "consensus"]<- "additive_Mla6_Bln1"
# bgh_spread_tables_Mla6_Bln1[is.na(bgh_spread_tables_Mla6_Bln1)]<- "additive_Mla6_Bln1"
bgh_spread_tables_Mla6_Bln1<- bgh_spread_tables_Mla6_Bln1[bgh_spread_tables_Mla6_Bln1$consensus!="bgh_additive_Mla6_Bln1",]

spread_tables_Mla6_Bln1<- rbind(spread_tables_Mla6_Bln1, bgh_spread_tables_Mla6_Bln1)

#pred_obs_epi<- pred_obs[pred_obs$gene %in% spread_tables_Mla6_Bln1$gene,]

gathered_tables_Mla6_Bln1<- gather(spread_tables_Mla6_Bln1, key = "time", value = "type", 2:8)
gathered_tables_Mla6_Bln1<- gathered_tables_Mla6_Bln1[!is.na(gathered_tables_Mla6_Bln1$type),]
#summary_barley_Mla6_Bln1<- gather(data = summary_barley_Mla6_Bln1, key = "time", value="freq",2:7)

pred_obs_epi<- merge(pred_obs, gathered_tables_Mla6_Bln1, by=c("gene", "time"))


pdf("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/global_epistasis/global_epistasis_epiGenes.pdf", width = 18, height = 12, fonts = "ArialMT", pointsize = 18)
palette="Set1"
ggplot(data = pred_obs_epi, aes_string(x = "pred_addi", y = "deviation", color = "time")) + theme_bw(base_size=40) +
  geom_point(size = 4) + coord_fixed()+ scale_color_brewer(palette = palette) +facet_wrap(~time)

ggplot(data = pred_obs_epi, aes_string(x = "pred_addi", y = "deviation", color = "time")) + theme_bw(base_size=40) +
  geom_point(size = 4) + coord_fixed()+ scale_color_brewer(palette = palette) +facet_wrap(~org) +
  geom_smooth(method = "lm", se = FALSE)

ggplot(data = pred_obs_epi[pred_obs_epi$org=="Barley",], aes_string(x = "pred_addi", y = "deviation", color = "time")) + theme_bw(base_size=40) +
  geom_point(size = 4) + coord_fixed()+ scale_color_brewer(palette = palette) 

ggplot(data = pred_obs_epi[pred_obs_epi$org=="Barley",], aes_string(x = "pred_addi", y = "deviation", color = "time")) + theme_bw(base_size=40) +
  geom_point(size = 4) + coord_fixed()+ scale_color_brewer(palette = palette) +facet_wrap(~time)

ggplot(data = pred_obs_epi[pred_obs_epi$org=="Bgh",], aes_string(x = "pred_addi", y = "deviation", color = "time")) + theme_bw(base_size=40) +
  geom_point(size = 4) + coord_fixed()+ scale_color_brewer(palette = palette)

ggplot(data = pred_obs_epi[pred_obs_epi$org=="Bgh",], aes_string(x = "pred_addi", y = "deviation", color = "time")) + theme_bw(base_size=40) +
  geom_point(size = 4) + coord_fixed()+ scale_color_brewer(palette = palette) +facet_wrap(~time)

dev.off()



global_epistasis_epi<- pred_obs_epi %>% group_by(org, time) %>% 
  summarise(slope=lm(deviation ~ pred_addi)$coefficients[2])
global_epistasis_epi$slope<- round(global_epistasis_epi$slope, digits = 4)

#model<- lm(pred_obs_deg$deviation ~ pred_obs_deg$pred_addi)
#model$coefficients[2]

global_epistasis_epi_all<- pred_obs_epi %>% group_by(org) %>% 
  summarise(slope=lm(deviation ~ pred_addi)$coefficients[2])
global_epistasis_epi_all$time<-"all"
global_epistasis_epi_all$slope<- round(global_epistasis_epi_all$slope, digits = 4)
global_epistasis_epi_all<-global_epistasis_epi_all[,c(1,3,2)]
global_epistasis_epi<- rbind(global_epistasis_epi, global_epistasis_epi_all)

write.csv(global_epistasis_epi, "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/global_epistasis/global_epsitasis_epigenes.csv")

