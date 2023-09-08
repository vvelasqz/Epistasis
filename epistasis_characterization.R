# description of epistasis genes for barley and bgh
library(tidyverse)

#barley
# spread_tables_Mla6_Bln1<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/spread_tables_Mla6_Bln1.csv", row.names = 1)
# spread_tables_Mla6_Rar3<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/spread_tables_Mla6_Rar3.csv", row.names = 1)
spread_tables_Mla6_Bln1<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/spread_tables_Mla6_Bln1.csv", row.names = 1)
spread_tables_Mla6_Rar3<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/spread_tables_Mla6_Rar3.csv", row.names = 1)
#spread_tables_Mla6_Bln1[spread_tables_Mla6_Bln1=="additive_Mla6_Bln1"]<-NA
#epistatic_spread_tables_Mla6_Bln1<- epistatic_spread_tables_Mla6_Bln1[rowSums(is.na(epistatic_spread_tables_Mla6_Bln1[,-1]))<6,]
#epistatic_spread_tables_Mla6_Bln1<- na.omit(epistatic_spread_tables_Mla6_Bln1[,-1])
# spread_tables_Mla6_Bln1$consensus<- apply(spread_tables_Mla6_Bln1[,2:7], 1, 
#                                           FUN=function(x){y=table(unlist(x)) 
#                                           paste(names(y)[which(y==max(y))], collapse = ";")})
# 
# spread_tables_Mla6_Bln1[spread_tables_Mla6_Bln1$consensus=="", "consensus"]<- "additive_Mla6_Bln1"
# spread_tables_Mla6_Bln1[is.na(spread_tables_Mla6_Bln1)]<- "additive_Mla6_Bln1"
spread_tables_Mla6_Bln1<- spread_tables_Mla6_Bln1[spread_tables_Mla6_Bln1$consensus!="", ]
consensus_barley_Mla6_Bln1<- data.frame(table(spread_tables_Mla6_Bln1$consensus))
consensus_barley_Mla6_Bln1<- consensus_barley_Mla6_Bln1 %>%  mutate(Type = strsplit(as.character(Var1), ";")) %>% unnest(Type) %>% 
  group_by(Type) %>% summarise(freq= sum(Freq))

# spread_tables_Mla6_Rar3$consensus<- apply(spread_tables_Mla6_Rar3[,2:7], 1, 
#                                           FUN=function(x){y=table(unlist(x)) 
#                                           paste(names(y)[which(y==max(y))], collapse = ";")})

consensus_barley_Mla6_Rar3<- data.frame(table(spread_tables_Mla6_Rar3$consensus))
consensus_barley_Mla6_Rar3<- consensus_barley_Mla6_Rar3 %>%  mutate(Type = strsplit(as.character(Var1), ";")) %>% unnest(Type) %>% 
  group_by(Type) %>% summarise(freq= sum(Freq))

consensus_barley_Mla6_Bln1$Type<- factor(consensus_barley_Mla6_Bln1$Type, 
                                         levels = c("additive_Mla6_Bln1", "symmetric", "masked", "suppression", "pseudo_masked", "positive", "negative"),
                                         labels = c("Additive", "Symmetric", "Masked", "Suppression", "Pseudo-masked", "Positive", "Negative"))
consensus_barley_Mla6_Rar3$Type<- factor(consensus_barley_Mla6_Rar3$Type, 
                                         levels = c( "unique_eff_Mla6", "unique_eff_Rar3", "equal_eff_Mla6_Rar3", "predominant_Mla6", "predominant_Rar3"),
                                         labels = c("Dominant Mla6", "Dominant Rar3", "Equal Mla6, Rar3", "Predominant Mla6", "Predominant Rar3"))

summary_barley_Mla6_Bln1<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/summary_barley_Mla6_Bln1.csv")
colnames(summary_barley_Mla6_Bln1)[1]<- "Type"
summary_barley_Mla6_Rar3<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/summary_barley_Mla6_Rar3.csv")
colnames(summary_barley_Mla6_Rar3)[1]<- "Type"

summary_barley_Mla6_Bln1<- gather(data = summary_barley_Mla6_Bln1, key = "time", value="freq",2:7)

summary_barley_Mla6_Rar3<- gather(data = summary_barley_Mla6_Rar3, key = "time", value="freq",2:7)
summary_barley_Mla6_Rar3[is.na(summary_barley_Mla6_Rar3$freq), "freq"]<- 0
summary_barley_Mla6_Rar3$Type<- factor(summary_barley_Mla6_Rar3$Type, 
                                       levels = c( "unique_eff_Mla6", "unique_eff_Rar3", "equal_eff_Mla6_Rar3", "predominant_Mla6", "predominant_Rar3"),
                                       labels = c("Dominant Mla6", "Dominant Rar3", "Equal Mla6, Rar3", "Predominant Mla6", "Predominant Rar3"))
summary_barley_Mla6_Bln1$Type<- factor(summary_barley_Mla6_Bln1$Type, 
                                       levels = c("additive_Mla6_Bln1", "symmetric", "masked", "suppression", "pseudo_masked", "positive", "negative"),
                                       labels = c("Additive", "Symmetric", "Masked", "Suppression", "Pseudo-masked", "Positive", "Negative"))

#for Bgh
bgh_spread_tables_Mla6_Bln1<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/bgh/HVR3/bgh_spread_tables_Mla6_Bln1.csv", row.names = 1)
bgh_spread_tables_Mla6_Rar3<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/bgh/HVR3/bgh_spread_tables_Mla6_Rar3.csv", row.names = 1)
# bgh_spread_tables_Mla6_Bln1[bgh_spread_tables_Mla6_Bln1=="additive_Mla6_Bln1"]<-NA
# bgh_spread_tables_Mla6_Bln1$consensus<- apply(bgh_spread_tables_Mla6_Bln1[,2:7], 1, 
#                                               FUN=function(x){y=table(unlist(x)) 
#                                               paste(names(y)[which(y==max(y))], collapse = ";")})
# 
# bgh_spread_tables_Mla6_Bln1[bgh_spread_tables_Mla6_Bln1$consensus=="", "consensus"]<- "additive_Mla6_Bln1"
# bgh_spread_tables_Mla6_Bln1[is.na(bgh_spread_tables_Mla6_Bln1)]<- "additive_Mla6_Bln1"

bgh_spread_tables_Mla6_Bln1<- bgh_spread_tables_Mla6_Bln1[bgh_spread_tables_Mla6_Bln1$consensus!="",]
consensus_bgh_Mla6_Bln1<- data.frame(table(bgh_spread_tables_Mla6_Bln1$consensus))
consensus_bgh_Mla6_Bln1<- consensus_bgh_Mla6_Bln1 %>%  mutate(Type = strsplit(as.character(Var1), ";")) %>% unnest(Type) %>% 
  group_by(Type) %>% summarise(freq= sum(Freq))

# bgh_spread_tables_Mla6_Rar3$consensus<- apply(bgh_spread_tables_Mla6_Rar3[,2:7], 1, 
#                                               FUN=function(x){y=table(unlist(x)) 
#                                               paste(names(y)[which(y==max(y))], collapse = ";")})

consensus_bgh_Mla6_Rar3<- data.frame(table(bgh_spread_tables_Mla6_Rar3$consensus))
consensus_bgh_Mla6_Rar3<- consensus_bgh_Mla6_Rar3 %>%  mutate(Type = strsplit(as.character(Var1), ";")) %>% unnest(Type) %>% 
  group_by(Type) %>% summarise(freq= sum(Freq))
consensus_bgh_Mla6_Bln1$Type<- factor(consensus_bgh_Mla6_Bln1$Type, 
                                      levels = c("bgh_additive_Mla6_Bln1", "bgh_symmetric", "bgh_masked", "bgh_suppression", "bgh_positive", "bgh_negative"),
                                      labels = c("Additive", "Symmetric", "Masked", "Suppression", "Positive", "Negative"))
consensus_bgh_Mla6_Rar3$Type<- factor(consensus_bgh_Mla6_Rar3$Type, 
                                      levels = c( "bgh_unique_eff_Mla6", "bgh_unique_eff_Rar3", "bgh_equal_eff_Mla6_Rar3", "bgh_predominant_Mla6"),
                                      labels = c("Dominant Mla6", "Dominant Rar3", "Equal Mla6, Rar3", "Predominant Mla6"))

summary_bgh_Mla6_Bln1<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/bgh/HVR3/summary_bgh_Mla6_Bln1.csv")
summary_bgh_Mla6_Rar3<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/bgh/HVR3/summary_bgh_Mla6_Rar3.csv")
colnames(summary_bgh_Mla6_Bln1)[1]<- "Type"
colnames(summary_bgh_Mla6_Rar3)[1]<- "Type"
summary_bgh_Mla6_Bln1<- gather(data = summary_bgh_Mla6_Bln1, key = "time", value="freq",2:7)
summary_bgh_Mla6_Rar3<- gather(data = summary_bgh_Mla6_Rar3, key = "time", value="freq",2:7)
summary_bgh_Mla6_Bln1$Type<- factor(summary_bgh_Mla6_Bln1$Type, 
                                    levels = c("bgh_additive_Mla6_Bln1", "bgh_symmetric", "bgh_masked", "bgh_suppression", "bgh_positive", "bgh_negative"),
                                    labels = c("Additive", "Symmetric", "Masked", "Suppression", "Positive", "Negative"))
summary_bgh_Mla6_Rar3[is.na(summary_bgh_Mla6_Rar3$freq), "freq"]<- 0
unique(summary_bgh_Mla6_Rar3$Type)
summary_bgh_Mla6_Rar3$Type<- factor(summary_bgh_Mla6_Rar3$Type, 
                                    levels = c( "bgh_unique_eff_Mla6", "bgh_unique_eff_Rar3", "bgh_equal_eff_Mla6_Rar3", "bgh_predominant_Mla6"),
                                    labels = c("Dominant Mla6", "Dominant Rar3", "Equal Mla6, Rar3", "Predominant Mla6"))

#barplot of patterns
#barley
#getPalette = c("#E69F00", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#56B4E9", "#000000")
#benchmarking<- read.csv("~/iowa_state/lab/Y2H publication/plos comp biol submission/second/figures_mla_version/benchmarking.csv")
#benchmarking$Reference<- factor(benchmarking$Reference, levels = c("Erffelinck et al. 2018", "Pashkova et al. 2016",   "Yachie et al. 2016", "Schlecht et al. 2017",   "Yang et al. 2018"))
pdf(paste0("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/barplots_R3.pdf"), width = 14, height = 12, fonts = "ArialMT", pointsize = 18)

getPalette  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(summary_barley_Mla6_Bln1, aes(x=Type, y=freq, group=time, fill=Type)) +
  geom_text(aes(label=freq),hjust=0.5, vjust=0, size=6) +geom_bar(stat="identity") + 
  scale_color_manual(values = getPalette, labels= c("Additive", "Symmetric", "Masked", "Suppression", "Pseudo-masked", "Positive", "Negative"))+ 
  scale_fill_manual(values = getPalette, labels= c("Additive", "Symmetric", "Masked", "Suppression", "Pseudo-masked", "Positive", "Negative"))+
  facet_wrap(~time, ncol = 3, scales = "free_y") + theme_grey(base_size = 30) + xlab("Type") +ylab("Number of barley genes") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.text.y = element_text(angle = 0), legend.position="bottom", 
        legend.text = element_text(margin = margin(r = 1.5, unit = "cm")),legend.margin = margin(l=0, unit="cm"))

ggplot(consensus_barley_Mla6_Bln1, aes(x=Type, y=freq, fill=Type)) +
  geom_text(aes(label=freq),hjust=0.5, vjust=0, size=15) +geom_bar(stat="identity") + 
  scale_color_manual(values = getPalette, labels= c("Additive", "Symmetric", "Masked", "Suppression", "Pseudo-masked", "Positive", "Negative"))+ 
  scale_fill_manual(values = getPalette, labels= c("Additive", "Symmetric", "Masked", "Suppression", "Pseudo-masked", "Positive", "Negative"))+
  theme_grey(base_size = 55) + xlab("Type") +ylab("Number of barley genes") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.text.y = element_text(angle = 0), legend.position="bottom", 
        legend.text = element_text(margin = margin(r = 1.5, unit = "cm")),legend.margin = margin(l=0, unit="cm"))

#Rar3
#getPalette = c("#FDE725FF","#21908CFF","#440154FF",  "#5DC863FF", "#3B528BFF")
getPalette = c("#00AFBB", "#E7B800", "#FC4E07", "#52854C", "#7570B3")
ggplot(summary_barley_Mla6_Rar3, aes(x=Type, y=freq, group=time, fill=Type)) +
  geom_text(aes(label=freq),hjust=0.5, vjust=0, size=6) +geom_bar(stat="identity") + 
  scale_color_manual(values = getPalette, labels= c("Dominant Mla6", "Dominant Rar3", "Equal Mla6, Rar3", "Predominant Mla6", "Predominant Rar3"))+ 
  scale_fill_manual(values = getPalette, labels= c("Dominant Mla6", "Dominant Rar3", "Equal Mla6, Rar3", "Predominant Mla6", "Predominant Rar3"))+
  facet_wrap(~time, ncol = 3, scales = "free_y") + theme_grey(base_size = 30) + xlab("Type") +ylab("Number of barley genes") + scale_y_continuous(trans='sqrt')+ 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.text.y = element_text(angle = 0), legend.position="bottom", 
        legend.text = element_text(margin = margin(r = 1.5, unit = "cm")), legend.margin = margin(l=0, unit="cm"))

ggplot(consensus_barley_Mla6_Rar3, aes(x=Type, y=freq, fill=Type)) +
  geom_text(aes(label=freq),hjust=0.5, vjust=0, size=15) +geom_bar(stat="identity") + 
  scale_color_manual(values = getPalette, labels= c("Dominant Mla6", "Dominant Rar3", "Equal Mla6, Rar3", "Predominant Mla6", "Predominant Rar3"))+ 
  scale_fill_manual(values = getPalette, labels= c("Dom Mla6", "Dom Rar3", "Equal", "Pred Mla6", "Pred Rar3"))+
  theme_grey(base_size = 55) + xlab("Type") +ylab("Number of barley genes") + scale_y_continuous(trans='sqrt')+ 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.text.y = element_text(angle = 0))


#bgh
getPalette  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")
names(getPalette)<- c("Additive", "Symmetric", "Masked", "Suppression", "Pseudo_masked", "Positive", "Negative")

ggplot(summary_bgh_Mla6_Bln1, aes(x=Type, y=freq, group=time, fill=Type)) +
  geom_text(aes(label=freq),hjust=0.5, vjust=0, size=6) +geom_bar(stat="identity") + 
  scale_color_manual(values = getPalette, labels= c("Additive", "Symmetric", "Masked", "Suppression",  "Positive", "Negative"))+ 
  scale_fill_manual(values = getPalette, labels= c("Additive", "Symmetric", "Masked", "Suppression",  "Positive", "Negative"))+
  facet_wrap(~time, ncol = 3, scales = "free_y") + theme_grey(base_size = 30) + xlab("Type") +ylab("Number of Bgh genes") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.text.y = element_text(angle = 0), legend.position="bottom", 
        legend.text = element_text(margin = margin(r = 1.5, unit = "cm")),legend.margin = margin(l=0, unit="cm"))

ggplot(consensus_bgh_Mla6_Bln1, aes(x=Type, y=freq, fill=Type)) +
  geom_text(aes(label=freq),hjust=0.5, vjust=0, size=15) +geom_bar(stat="identity") + 
  scale_color_manual(values = getPalette, labels= c("Additive", "Symmetric", "Masked", "Suppression", "Positive", "Negative"))+ 
  scale_fill_manual(values = getPalette, labels= c("Additive", "Symmetric", "Masked", "Suppression", "Positive", "Negative"))+
  theme_grey(base_size = 55) + xlab("Type") +ylab("Number of Bgh genes") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.text.y = element_text(angle = 0), legend.position="bottom", 
        legend.text = element_text(margin = margin(r = 1.5, unit = "cm")),legend.margin = margin(l=0, unit="cm"))

#Rar3
#getPalette = c("#FDE725FF","#21908CFF","#440154FF",  "#5DC863FF", "#3B528BFF")
getPalette = c("#00AFBB", "#E7B800", "#FC4E07", "#52854C", "#7570B3")
ggplot(summary_bgh_Mla6_Rar3, aes(x=Type, y=freq, group=time, fill=Type)) +
  geom_text(aes(label=freq),hjust=0.5, vjust=0, size=6) +geom_bar(stat="identity") + 
  scale_color_manual(values = getPalette, labels= c("Dominant Mla6", "Dominant Rar3", "Equal Mla6, Rar3", "Predominant Mla6"))+ 
  scale_fill_manual(values = getPalette, labels= c("Dominant Mla6", "Dominant Rar3", "Equal Mla6, Rar3", "Predominant Mla6"))+
  facet_wrap(~time, ncol = 3, scales = "free_y") + theme_grey(base_size = 30) + xlab("Type") +ylab("Number of Bgh genes") + scale_y_continuous(trans='sqrt')+ 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.text.y = element_text(angle = 0), legend.position="bottom", 
        legend.text = element_text(margin = margin(r = 1.5, unit = "cm")),legend.margin = margin(l=0, unit="cm"))

ggplot(consensus_bgh_Mla6_Rar3, aes(x=Type, y=freq, fill=Type)) +
  geom_text(aes(label=freq),hjust=0.5, vjust=0, size=15) +geom_bar(stat="identity") + 
  scale_color_manual(values = getPalette, labels= c("Dominant Mla6", "Dominant Rar3", "Equal Mla6, Rar3", "Predominant Mla6"))+ 
  scale_fill_manual(values = getPalette, labels= c("Dom Mla6", "Dom Rar3", "Equal", "Pred Mla6"))+
  theme_grey(base_size = 55) + xlab("Type") +ylab("Number of Bgh genes") + scale_y_continuous(trans='sqrt')+ 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.text.y = element_text(angle = 0))

dev.off()




#create the lists for chr location
#directory<- "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/chr analysis/R2/Genelists for positioning/"
directory<- "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/chr analysis/R3/Genelists for positioning/"
#annotation<- read.table("~/iowa_state/lab/RNAseq/DESeq2 HvR2/annotation_HvR2.txt", sep=';', stringsAsFactors = F, quote = "\"")[,1:2]
annotation<- read.csv("~/iowa_state/lab/genome annotation files/tritex R3/HR3_annotation.csv")[,c(3,6)]
colnames(annotation)<-c("gene","description")
epistatic_spread_tables_Mla6_Bln1<- spread_tables_Mla6_Bln1[spread_tables_Mla6_Bln1$consensus!="additive_Mla6_Bln1", ]
epistatic_spread_tables_Mla6_Bln1<- merge(epistatic_spread_tables_Mla6_Bln1, annotation, by="gene", all.x=T)
for(type in unique(epistatic_spread_tables_Mla6_Bln1$consensus)){
  table<- epistatic_spread_tables_Mla6_Bln1[epistatic_spread_tables_Mla6_Bln1$consensus==type, c(1,9,2:8)]
  write.csv(table, paste0(directory, "Barley_Mla6_Bln1_", type, "_Genelist.csv"))
}

spread_tables_Mla6_Rar3<- merge(spread_tables_Mla6_Rar3, annotation, by="gene", all.x=T)
for(type in unique(spread_tables_Mla6_Rar3$consensus)){
  table<- spread_tables_Mla6_Rar3[spread_tables_Mla6_Rar3$consensus==type, c(1,9,2:8)]
  write.csv(table, paste0(directory, "Barley_Mla6_Rar3_", type, "_Genelist.csv"))
}
write.csv(c(paste0("Barley_Mla6_Bln1_",unique(epistatic_spread_tables_Mla6_Bln1$consensus)),
            paste0("Barley_Mla6_Rar3_",unique(spread_tables_Mla6_Rar3$consensus))),
          "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/chr analysis/R3/trts_barley.csv")

# write.csv(c(paste0("Barley_Mla6_Bln1_",unique(epistatic_spread_tables_Mla6_Bln1$consensus)),
#             paste0("Barley_Mla6_Rar3_",unique(spread_tables_Mla6_Rar3$consensus))),
#           "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/chr analysis/R2/trts_barley.csv")




# spread_tables_Mla6_Bln1$consensus<- apply(spread_tables_Mla6_Bln1[,2:7], 1, 
#                                           FUN=function(x){y=table(unlist(x)) 
#                                           paste(names(y)[which(y==max(y))], collapse = ";")})

# y= table(unlist(spread_tables_Mla6_Bln1[2,2:7]))
# y
# names(y)[which(y==max(y))]



# #position analysis
# #folder<- "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/chr analysis/R2/"
# folder<- "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/chr analysis/R3/"
# #start by creating the position file from gff 
# barley_gff <- read.table(file="~/iowa_state/lab/genome annotation files/tritex R3/Hv_Morex.pgsb.Jul2020.gff3",
#                          sep = '\t',header = F, quote = "\"", fill=TRUE, stringsAsFactors = F, skip = 9)
# 
# # barley_gff <- read.table(file="~/iowa_state/lab/genome annotation files/tritex R2/Barley_Morex_V2_gene_annotation_PGSB.HC.gff3", 
# #                       sep = '\t',header = F, quote = "\"", fill=TRUE, stringsAsFactors = F, skip = 9)
# 
# barley_gff<- barley_gff[barley_gff$V3=="gene",c(1,4,5,7,9)]
# barley_gff$V9<- sapply(X = barley_gff$V9, FUN=function(x) substr(x, 4,28))
# colnames(barley_gff)<- paste0("V", 1:5)
# 
# write.csv(barley_gff,paste0(folder, "Horvu_table.csv"))
# 
# horvu_V2 <- read.table('~/iowa_state/lab/genome annotation files/mloc to horvu/Horvu_table.tsv',
#                     stringsAsFactors = F)
# #plot positions
# #http://bioconductor.org/packages/devel/bioc/vignettes/profileplyr/inst/doc/profileplyr.html
# #https://compgenomr.github.io/book/visualizing-and-summarizing-genomic-intervals.html#visualizing-intervals-on-a-locus-of-interest
# 
# library(tidyverse)
# library(GenomicRanges)
# horvu <- read.csv(paste0(folder, "Horvu_table.csv"),row.names = 1, stringsAsFactors = F)
# len<- horvu %>% group_by(V1) %>% summarise(lenght= max(V2, V3))
# names(horvu)<- (c("chr","start","end","strand","seqname"))
# Seqinfo(as.character(len$V1), seqlengths=len$lenght, isCircular=NA, genome=NA)
# horvu_ranges<- makeGRangesFromDataFrame(horvu,keep.extra.columns=TRUE,ignore.strand=FALSE,
#                                         seqinfo=Seqinfo(as.character(len$V1), seqlengths=len$lenght, isCircular=NA, genome=NA),
#                                         seqnames.field="chr",start.field="start",end.field="end",
#                                         strand.field="strand",starts.in.df.are.0based=FALSE)
# 
# 
# 
# bins<- data.frame(chromosome= rep(unlist(strsplit(len$V1, "r"))[-which(unlist(strsplit(len$V1, "r"))=="ch")], ceiling(len$lenght/1e6)), 
#                   chr= rep(len$V1, ceiling(len$lenght/1e6)),
#                   number= unlist(sapply(len$lenght, function(x) sprintf("%03d",0:floor(x/1e6)))),
#                   start= unlist(sapply(len$lenght, function(x) 0:floor(x/1e6)*1e6+1)), 
#                   end= unlist(sapply(len$lenght, function(x) 0:floor(x/1e6)*1e6+1e6)),
#                   strand="+")
# 
# bins$number<- as.character(bins$number)
# bins$bin_1Mb<- sapply(sapply(bins$number, function(x) strsplit(x, "")), function(x) paste(x, collapse = ".")) 
# bins$bin_1Mb<- apply(bins, 1, function(x) paste0(x[1],".", x[7])) 
# bins[sapply(1:nrow(len), function(x) sum(ceiling(len$lenght/1e6)[1:x])), "end"]<- len$lenght
# 
# bins_ranges<- makeGRangesFromDataFrame(bins,keep.extra.columns=TRUE,ignore.strand=FALSE,
#                                        seqinfo=Seqinfo(as.character(len$V1), seqlengths=len$lenght, isCircular=NA, genome=NA),
#                                        seqnames.field="chr",start.field="start",end.field="end",
#                                        strand.field="strand",starts.in.df.are.0based=FALSE)
# 
# complete_overlap = findOverlaps(query = horvu_ranges, subject = bins_ranges, type = "any", ignore.strand=T)
# complete_overlap_df = data.frame(horvu[queryHits(complete_overlap),], bins[subjectHits(complete_overlap),])
# complete_overlap_df$bin_10Mb<- gsub('.{0,2}$', '', complete_overlap_df$bin_1Mb) 
# complete_overlap_df$bin_100Mb<- gsub('.{0,4}$', '', complete_overlap_df$bin_1Mb) 
# #colnames(complete_overlap_df)[11:13]<- c("bin_1Mb", "bin_10Mb", "bin_100Mb")
# #write.csv(complete_overlap_df, paste0(folder, "Binned_genes_HR2.csv"), row.names = F)
# write.csv(complete_overlap_df, paste0(folder, "Binned_genes_HR3.csv"), row.names = F)
# 
# genes_bin<-data.frame(table(complete_overlap_df$bin_1Mb))
# bins<- merge(bins, genes_bin, by.x="bin_1Mb", by.y="Var1", all.x=T)
# bins[is.na(bins$Freq),"Freq"]<-0
# bins$bin_10Mb<- gsub('.{0,2}$', '', bins$bin_1Mb) 
# bins$bin_100Mb<- gsub('.{0,4}$', '', bins$bin_1Mb) 
# 
# bins_1Mb<- bins %>% group_by(bin_1Mb) %>% summarise(chromosome=unique(chromosome), start=(min(start)), end=max(end), Freq=sum(Freq))
# bins_10Mb<- bins %>% group_by(bin_10Mb) %>% summarise(chromosome=unique(chromosome), start=(min(start)), end=max(end), Freq=sum(Freq))
# bins_100Mb<- bins %>% group_by(bin_100Mb) %>% summarise(chromosome=unique(chromosome), start=(min(start)), end=max(end), Freq=sum(Freq))
# bins_chr<- bins %>% group_by(chromosome) %>% summarise(start=(min(start)), end=max(end), Freq=sum(Freq))
# 
# write.csv(bins_1Mb, paste0(folder, "bins_1Mb.csv"), row.names = F)
# write.csv(bins_10Mb, paste0(folder, "bins_10Mb.csv"), row.names = F)
# write.csv(bins_100Mb, paste0(folder, "bins_100Mb.csv"), row.names = F)
# write.csv(bins_chr, paste0(folder, "bins_chrs.csv"), row.names = F)
# 
# #HG test
# HG_test<- function(trts, bin_table, overlap_table){
#   #table trts by bins
#   hypergeom_test <- data.frame(matrix(ncol = 1+2*length(trts), nrow = nrow(bin_table))) 
#   hypergeom_test[,1]<- bin_table[,1]
#   colnames(hypergeom_test)<- c("binID", unlist(sapply(trts, function(x) c(paste0(x,"_pvalue"), paste0(x,"_padj")))))
#   
#   for(m in 1:length(trts)){
#     epi_gene_list<- read.csv(paste0(folder,"Genelists for positioning/",trts[m], "_Genelist.csv"), 
#                             stringsAsFactors = F, row.names = 1)[,"gene"]
#     overlap_table_filtered<- overlap_table[overlap_table$seqname %in% epi_gene_list, c("seqname", colnames(bin_table)[1])]
#     DE_genes_bins<-data.frame(table(overlap_table_filtered[,colnames(bin_table)[1]]))
#     overlap_table_filtered<- merge(overlap_table_filtered, DE_genes_bins, by.x=colnames(bin_table)[1], by.y="Var1", all.x=T)
#     DE_bins<- bin_table[unlist(bin_table[,1]) %in% DE_genes_bins$Var1,]
#     
#     #total DE in bin, total DE trt, sum not DE in bins listed in trt, genes in bin
#     hypergeom_test[match(DE_genes_bins$Var1, unlist(hypergeom_test$binID)), 
#                    paste0(trts[m],"_pvalue")]<- phyper(DE_genes_bins$Freq, length(epi_gene_list), 
#                                                        sum(DE_bins$Freq - DE_genes_bins$Freq), 
#                                                        DE_bins$Freq, lower.tail = FALSE)
#     
#     hypergeom_test[unlist(hypergeom_test$binID) %in% DE_genes_bins$Var1,paste0(trts[m],"_padj")]<- p.adjust(hypergeom_test[unlist(hypergeom_test$binID) %in% DE_genes_bins$Var1,paste0(trts[m],"_pvalue")], method = "BH")
#   }
#   hypergeom_test<- hypergeom_test[rowSums(is.na(hypergeom_test[,-1])) != ncol(hypergeom_test[,-1]), ]
#   return(hypergeom_test)
# }
# 
# #test for all bin tables
# chrs = list.files(path = paste0(folder,"Genelists for positioning/"), pattern = "csv")
# #trts<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/chr analysis/R2/trts_barley.csv", row.names = 1)
# trts<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/chr analysis/R3/trts_barley.csv", row.names = 1)
# trts<- trts$x
# #trts<- gsub('.{0,13}$', '', chrs) 
# #overlap_table<- read.csv(paste0(folder, "Binned_genes.csv")
# #overlap_table<- read.csv(paste0(folder, "Binned_genes_HR3.csv")
# 
# overlap_table<- complete_overlap_df
# 
# # bins_1Mb<- read.csv( paste0(folder, "bins_1Mb.csv")
# # bins_10Mb<-read.csv(paste0(folder, "bins_10Mb.csv")
# # bins_100Mb<-read.csv( paste0(folder, "bins_100Mb.csv")
# # bins_chr<-read.csv( paste0(folder, "bins_chrs.csv")
# 
# 
# HG_1Mb<- HG_test(trts, bins_1Mb, overlap_table)
# write.csv(HG_1Mb, paste0(folder, "hypergeometric_test_bins_1Mb.csv"), row.names = F)
# 
# HG_10Mb<- HG_test(trts, bins_10Mb, overlap_table)
# write.csv(HG_10Mb, paste0(folder, "hypergeometric_test_bins_10Mb.csv"))
# 
# HG_100Mb<- HG_test(trts, bins_100Mb, overlap_table)
# write.csv(HG_100Mb, paste0(folder, "hypergeometric_test_bins_100Mb.csv"))
# 
# HG_chrs<- HG_test(trts, bins_chr, overlap_table)
# write.csv(HG_chrs, paste0(folder, "hypergeometric_test_bins_chrs.csv"))
# 
# 
# #Make the plots
# #https://www.bioconductor.org/packages/release/bioc/vignettes/ggbio/inst/doc/ggbio.pdf
# #BiocManager::install("ggbio")
# 
# library(ggbio)
# #Now a loop
# p <- autoplot(seqinfo(horvu_ranges), layout = "karyogram")
# #pdf(paste0(folder, "Barley_Mla6_Bln1_chr_position_by_group.pdf"), width = 10, height = 10, fonts = "ArialMT", pointsize = 12)
# pdf(paste0(folder, "Barley_chr_position_by_group.pdf"), width = 10, height = 10, fonts = "ArialMT", pointsize = 12)
# for(m in 1:length(chrs)){
#   name<- gsub('.{0,13}$', '', chrs[m])
#   coords<- read.csv(paste0(folder,"Genelists for positioning/",chrs[m]), stringsAsFactors = F, row.names = 1)
#   coords<- merge(horvu, coords, by.x="seqname",by.y="gene")
#   coords_ann<-makeGRangesFromDataFrame(coords, keep.extra.columns=T, seqnames.field= "chr",start.field="start",
#                                        seqinfo=Seqinfo(as.character(len$V1)), end.field="end", strand.field="strand")
#   
#   print(p + layout_karyogram(coords_ann, color="red") +labs(title = name))
#   # if(!grepl("dmDiff", name)){
#   #   print(p + layout_karyogram(coords_ann[coords_ann$Direction=="Pos",], color="red") +
#   #           layout_karyogram(coords_ann[coords_ann$Direction=="Neg",], color="blue")+labs(title = name))
#   # }else{
#   #   print(p + layout_karyogram(coords_ann, color="purple") +labs(title = name))
#   # }
# }
# dev.off()
# 
# 
# #Fig 4a position
# #define the significant hotspots and plot coloring by pattern
# folder<- "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/chr analysis/R3/"
# # HG_1Mb<- read.csv(paste0(folder, "hypergeometric_test_bins_1Mb.csv"))
# # HG_10Mb<- read.csv(paste0(folder, "hypergeometric_test_bins_10Mb.csv"), row.names = 1)
# HG_100Mb<- read.csv(paste0(folder, "hypergeometric_test_bins_100Mb.csv"), row.names = 1)
# #HG_chrs<- read.csv(paste0(folder, "hypergeometric_test_bins_chrs.csv"), row.names = 1)
# 
# # type=c("additive_Mla6_Bln1", "symmetric", "masked", "suppression", "pseudo_masked", "positive", "negative",
# #        "unique_eff_Mla6", "unique_eff_Rar3", "equal_eff_Mla6_Rar3", "predominant_Mla6", "predominant_Rar3")
# type=c("additive_Mla6_Bln1", "symmetric", "masked", "suppression", "pseudo_masked", "positive", "negative")
# trts<- paste0("Barley_Mla6_Bln1_",type, "_padj")
# 
# HG_100Mb_mla6_bln1<- HG_100Mb[, c(1, grep(paste(trts, collapse = "|"), colnames(HG_100Mb)))]
# #thr 0.005
# HG_100Mb_mla6_bln1<- HG_100Mb_mla6_bln1[rowSums(HG_100Mb_mla6_bln1[,-1] <0.005, na.rm = T)>0,
#                                         c(T,colSums(HG_100Mb_mla6_bln1[,-1] <0.005, na.rm = T)>0)]
# sig_trts<- sapply(colnames(HG_100Mb_mla6_bln1)[-1], FUN=function(x)strsplit(x,"_")[[1]][4])
# colnames(HG_100Mb_mla6_bln1)[-1]<-sig_trts
# HG_100Mb_mla6_bln1<- melt(HG_100Mb_mla6_bln1, 1, na.rm = T)
# HG_100Mb_mla6_bln1<- HG_100Mb_mla6_bln1[HG_100Mb_mla6_bln1$value< 0.005,]
# 
# binned_genes<-read.csv(paste0(folder, "Binned_genes_HR3.csv"))
# spread_tables_Mla6_Bln1<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/spread_tables_Mla6_Bln1.csv", row.names = 1)
# spread_tables_Mla6_Bln1<- spread_tables_Mla6_Bln1[spread_tables_Mla6_Bln1$consensus %in% sig_trts, c("gene","consensus")]
# spread_tables_Mla6_Bln1<- merge(spread_tables_Mla6_Bln1, binned_genes[,c(5,14)], by.x="gene", by.y="seqname", all.x=T)
# sig_genes_Mla6_Bln1<- merge(spread_tables_Mla6_Bln1, HG_100Mb_mla6_bln1, by.x=c("consensus","bin_100Mb"), by.y=c("variable", "binID"))
# 
# write.csv(sig_genes_Mla6_Bln1, paste0(folder, "Barley_chr_position_significant.csv"))
# 
# library(tidyverse)
# library(GenomicRanges)
# horvu <- read.csv(paste0(folder, "Horvu_table.csv"),row.names = 1, stringsAsFactors = F)
# len<- horvu %>% group_by(V1) %>% summarise(lenght= max(V2, V3))
# names(horvu)<- (c("chr","start","end","strand","seqname"))
# Seqinfo(as.character(len$V1), seqlengths=len$lenght, isCircular=NA, genome=NA)
# horvu_ranges<- makeGRangesFromDataFrame(horvu,keep.extra.columns=TRUE,ignore.strand=FALSE,
#                                         seqinfo=Seqinfo(as.character(len$V1), seqlengths=len$lenght, isCircular=NA, genome=NA),
#                                         seqnames.field="chr",start.field="start",end.field="end",
#                                         strand.field="strand",starts.in.df.are.0based=FALSE)
# 
# #Make the plots
# #https://www.bioconductor.org/packages/release/bioc/vignettes/ggbio/inst/doc/ggbio.pdf
# #BiocManager::install("ggbio")
# 
# library(ggbio)
# #Now a loop
# p <- autoplot(seqinfo(horvu_ranges), layout = "karyogram")
# chrs<- sig_trts
# colors  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# names(colors)<- c("additive_Mla6_Bln1", "symmetric", "masked", "suppression", "pseudo_masked", "positive", "negative")
# 
# for(m in 1:length(chrs)){
#   name<- chrs[m]
#   coords<- sig_genes_Mla6_Bln1[sig_genes_Mla6_Bln1$consensus==chrs[m],]
#   coords<- merge(horvu, coords, by.x="seqname",by.y="gene")
#   coords_ann<-makeGRangesFromDataFrame(coords, keep.extra.columns=T, seqnames.field= "chr",start.field="start",
#                                        seqinfo=Seqinfo(as.character(len$V1)), end.field="end", strand.field="strand")
#   p<-p + layout_karyogram(coords_ann, color=colors[chrs[m]])
# }
# 
# pdf(paste0(folder, "Barley_chr_position_significant.pdf"), width = 5, height = 3, fonts = "ArialMT", pointsize = 52)
# print(p)
# dev.off()
# 
# 
# #mla6 rar3
# folder<- "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/chr analysis/R3/"
# HG_100Mb<- read.csv(paste0(folder, "hypergeometric_test_bins_100Mb.csv"), row.names = 1)
# type=c("unique_eff_Mla6", "unique_eff_Rar3", "equal_eff_Mla6_Rar3", "predominant_Mla6", "predominant_Rar3")
# trts<- paste0("Barley_Mla6_Rar3_",type, "_padj")
# 
# HG_100Mb_mla6_rar3<- HG_100Mb[, c(1, grep(paste(trts, collapse = "|"), colnames(HG_100Mb)))]
# #thr 0.005
# HG_100Mb_mla6_rar3<- HG_100Mb_mla6_rar3[rowSums(HG_100Mb_mla6_rar3[,-1] <0.005, na.rm = T)>0,
#                                         c(T,colSums(HG_100Mb_mla6_rar3[,-1] <0.005, na.rm = T)>0)]
# sig_trts<- c("unique_eff_Mla6", "equal_eff_Mla6_Rar3", "predominant_Mla6")
# colnames(HG_100Mb_mla6_rar3)[-1]<-sig_trts
# HG_100Mb_mla6_rar3<- melt(HG_100Mb_mla6_rar3, 1, na.rm = T)
# HG_100Mb_mla6_rar3<- HG_100Mb_mla6_rar3[HG_100Mb_mla6_rar3$value< 0.005,]
# 
# binned_genes<-read.csv(paste0(folder, "Binned_genes_HR3.csv"))
# spread_tables_Mla6_Rar3<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/spread_tables_Mla6_Rar3.csv", row.names = 1)
# spread_tables_Mla6_Rar3<- spread_tables_Mla6_Rar3[spread_tables_Mla6_Rar3$consensus %in% sig_trts, c("gene","consensus")]
# spread_tables_Mla6_Rar3<- merge(spread_tables_Mla6_Rar3, binned_genes[,c(5,14)], by.x="gene", by.y="seqname", all.x=T)
# sig_genes_Mla6_Rar3<- merge(spread_tables_Mla6_Rar3, HG_100Mb_mla6_rar3, by.x=c("consensus","bin_100Mb"), by.y=c("variable", "binID"))
# 
# write.csv(sig_genes_Mla6_Rar3, paste0(folder, "Barley_chr_position_significant_rar3.csv"))
# 
# library(tidyverse)
# library(GenomicRanges)
# horvu <- read.csv(paste0(folder, "Horvu_table.csv"),row.names = 1, stringsAsFactors = F)
# len<- horvu %>% group_by(V1) %>% summarise(lenght= max(V2, V3))
# names(horvu)<- (c("chr","start","end","strand","seqname"))
# Seqinfo(as.character(len$V1), seqlengths=len$lenght, isCircular=NA, genome=NA)
# horvu_ranges<- makeGRangesFromDataFrame(horvu,keep.extra.columns=TRUE,ignore.strand=FALSE,
#                                         seqinfo=Seqinfo(as.character(len$V1), seqlengths=len$lenght, isCircular=NA, genome=NA),
#                                         seqnames.field="chr",start.field="start",end.field="end",
#                                         strand.field="strand",starts.in.df.are.0based=FALSE)
# 
# library(ggbio)
# #Now a loop
# p <- autoplot(seqinfo(horvu_ranges), layout = "karyogram")
# chrs<- sig_trts[-3]
# colors = c("#00AFBB", "#E7B800", "#FC4E07", "#52854C", "#7570B3")
# names(colors)<- c("unique_eff_Mla6", "unique_eff_Rar3", "equal_eff_Mla6_Rar3", "predominant_Mla6", "predominant_Rar3")
# 
# for(m in 1:length(chrs)){
#   name<- chrs[m]
#   coords<- sig_genes_Mla6_Rar3[sig_genes_Mla6_Rar3$consensus==chrs[m],]
#   coords<- merge(horvu, coords, by.x="seqname",by.y="gene")
#   coords_ann<-makeGRangesFromDataFrame(coords, keep.extra.columns=T, seqnames.field= "chr",start.field="start",
#                                        seqinfo=Seqinfo(as.character(len$V1)), end.field="end", strand.field="strand")
#   p<-p + layout_karyogram(coords_ann, color=colors[chrs[m]])
# }
# 
# pdf(paste0(folder, "Barley_chr_position_significant_mla6_rar3.pdf"), width = 5, height = 3, fonts = "ArialMT", pointsize = 52)
# print(p)
# dev.off()
# 
# 
# 




#go analysis

#Fig 3b
#Go analysis ############################
#barley
spread_tables_Mla6_Bln1<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/spread_tables_Mla6_Bln1.csv", row.names = 1)
spread_tables_Mla6_Rar3<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/spread_tables_Mla6_Rar3.csv", row.names = 1)
# spread_tables_Mla6_Bln1[spread_tables_Mla6_Bln1=="additive_Mla6_Bln1"]<-NA
# spread_tables_Mla6_Bln1$consensus<- apply(spread_tables_Mla6_Bln1[,2:7], 1, 
#                                           FUN=function(x){y=table(unlist(x)) 
#                                           paste(names(y)[which(y==max(y))], collapse = ";")})
# spread_tables_Mla6_Bln1[spread_tables_Mla6_Bln1$consensus=="", "consensus"]<- "additive_Mla6_Bln1"
# spread_tables_Mla6_Bln1[is.na(spread_tables_Mla6_Bln1)]<- "additive_Mla6_Bln1"
# 
# spread_tables_Mla6_Rar3$consensus<- apply(spread_tables_Mla6_Rar3[,2:7], 1, 
#                                           FUN=function(x){y=table(unlist(x)) 
#                                           paste(names(y)[which(y==max(y))], collapse = ";")})

spread_tables_barley <- rbind(spread_tables_Mla6_Bln1, spread_tables_Mla6_Rar3)

library(clusterProfiler)
#GO_terms<- read.delim("~/iowa_state/lab/GRN_reconstruction/HORVUR2 version analysis/ensemble of GRN/GO analyses/data/annotation_HvR2.txt", sep = ";", stringsAsFactors = F, header = F, fill=T)
GO_terms<- read.csv("~/iowa_state/lab/genome annotation files/tritex R3/HR3_annotation.csv", row.names = 1)[,c(2,5,6)]
colnames(GO_terms)<- c("gene", "description", "go_terms")
s <- strsplit(GO_terms$go_terms, split = ",")
GO_terms<- data.frame(gene = rep(GO_terms$gene, sapply(s, length)), description = rep(GO_terms$description, sapply(s, length)),GO = unlist(s))
rm(s)
term2gene=GO_terms[, c("GO", "gene")]
term2name=GO_terms[, c("GO", "description")]
#uniqueGO<- read.csv("~/iowa_state/lab/GRN_reconstruction/HORVUR2 version analysis/ensemble of GRN/GO analyses/uniqueGO.csv")
uniqueGO<- read.csv("~/iowa_state/lab/genome annotation files/tritex R3/HR3_uniqueGO.csv", row.names = 1)

#now DE targets per timepoint all fc
unique(spread_tables_barley$consensus)
type=c("additive_Mla6_Bln1", "symmetric", "masked", "suppression", "pseudo_masked", "positive", "negative",
       "unique_eff_Mla6", "unique_eff_Rar3", "equal_eff_Mla6_Rar3", "predominant_Mla6", "predominant_Rar3")
type %in% unique(spread_tables_barley$consensus)

go_result<- data.frame()

for (t in type){
  DE_gene_list<- spread_tables_barley[spread_tables_barley$consensus==t,"gene"]
  x_horvu <- enricher(DE_gene_list, TERM2GENE=term2gene, TERM2NAME=term2name)
  if(!is.null(x_horvu)){
    if(nrow(x_horvu)>0){
      x_horvu<- merge(as.data.frame(x_horvu), uniqueGO, by.x="ID", by.y="GO")
      x_horvu$group<- t
      go_result<- rbind(go_result, x_horvu)
    }}
  }

go_result<- go_result[!duplicated(go_result),]
write.csv(go_result, "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/GO/R3/GO_epistasis_barley.csv")

length(unique(go_result$ID))
length(unique(go_result$Description))

#go analysis per timepoint 

go_result_time<- data.frame()

for (time in c("t0","t16","t20","t24","t32","t48")){
  for (t in type){
    if(t %in% spread_tables_barley[,time]){
      DE_gene_list<- spread_tables_barley[spread_tables_barley[,time]==t,"gene"]
      #DE_gene_list<- spread_tables_barley[spread_tables_barley$consensus==t,"gene"]
      x_horvu <- enricher(DE_gene_list, TERM2GENE=term2gene, TERM2NAME=term2name)
      if(!is.null(x_horvu)){
        if(nrow(x_horvu)>0){
          x_horvu<- merge(as.data.frame(x_horvu), uniqueGO, by.x="ID", by.y="GO")
          x_horvu$type<- t
          x_horvu$group<- paste0(time, "_", t)
          x_horvu$time<- time
          go_result_time<- rbind(go_result_time, x_horvu)
        }
      }
    }
  }
}


go_result_time<- go_result_time[!duplicated(go_result_time),]
write.csv(go_result_time, "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/GO/R3/GO_epistasis_barley_timepoints.csv")

length(unique(go_result_time$ID))
length(unique(go_result_time$Description))




#GO Bgh
bgh_spread_tables_Mla6_Bln1<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/bgh/HVR3/bgh_spread_tables_Mla6_Bln1.csv", row.names = 1)
bgh_spread_tables_Mla6_Rar3<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/bgh/HVR3/bgh_spread_tables_Mla6_Rar3.csv", row.names = 1)
# bgh_spread_tables_Mla6_Bln1[bgh_spread_tables_Mla6_Bln1=="additive_Mla6_Bln1"]<-NA
# bgh_spread_tables_Mla6_Bln1$consensus<- apply(bgh_spread_tables_Mla6_Bln1[,2:7], 1, 
#                                           FUN=function(x){y=table(unlist(x)) 
#                                           paste(names(y)[which(y==max(y))], collapse = ";")})
# bgh_spread_tables_Mla6_Bln1[bgh_spread_tables_Mla6_Bln1$consensus=="", "consensus"]<- "additive_Mla6_Bln1"
# bgh_spread_tables_Mla6_Bln1[is.na(bgh_spread_tables_Mla6_Bln1)]<- "additive_Mla6_Bln1"
# 
# bgh_spread_tables_Mla6_Rar3$consensus<- apply(bgh_spread_tables_Mla6_Rar3[,2:7], 1, 
#                                           FUN=function(x){y=table(unlist(x)) 
#                                           paste(names(y)[which(y==max(y))], collapse = ";")})

bgh_spread_tables <- rbind(bgh_spread_tables_Mla6_Bln1, bgh_spread_tables_Mla6_Rar3)

library(clusterProfiler)
#GO_terms<- read.delim("~/iowa_state/lab/GRN_reconstruction/HORVUR2 version analysis/ensemble of GRN/GO analyses/data/annotation_HvR2.txt", sep = ";", stringsAsFactors = F, header = F, fill=T)
GO_terms<- read.csv("~/iowa_state/lab/genome annotation files/ef2 blumeria genome/dh14_GO_terms.csv", row.names = 1)
term2gene=GO_terms[, c("GO", "gene")]
term2name=GO_terms[, c("GO", "description")]
#uniqueGO<- read.csv("~/iowa_state/lab/GRN_reconstruction/HORVUR2 version analysis/ensemble of GRN/GO analyses/uniqueGO.csv")
uniqueGO<- read.csv("~/iowa_state/lab/genome annotation files/ef2 blumeria genome/dh14_uniqueGO.csv", row.names = 1)

unique(bgh_spread_tables$consensus)
type=c("bgh_additive_Mla6_Bln1", "bgh_symmetric", "bgh_masked", "bgh_suppression", "bgh_pseudo_masked", "bgh_positive", "bgh_negative",
       "bgh_unique_eff_Mla6", "bgh_unique_eff_Rar3", "bgh_equal_eff_Mla6_Rar3", "bgh_predominant_Mla6", "bgh_predominant_Rar3")
type %in% unique(bgh_spread_tables$consensus)

go_result<- data.frame()

for (t in type){
  if(t %in% bgh_spread_tables$consensus){
    DE_gene_list<- bgh_spread_tables[bgh_spread_tables$consensus==t,"gene"]
    x_horvu <- enricher(DE_gene_list, TERM2GENE=term2gene, TERM2NAME=term2name)
    if(!is.null(x_horvu)){
      if(nrow(x_horvu)>0){
        x_horvu<- merge(as.data.frame(x_horvu), uniqueGO, by.x="ID", by.y="GO")
        x_horvu$group<- t
        go_result<- rbind(go_result, x_horvu)
      }
    }
  }
}

go_result<- go_result[!duplicated(go_result),]
write.csv(go_result, "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/GO/R3/GO_epistasis_bgh.csv")

length(unique(go_result$ID))
length(unique(go_result$Description))

#interproscan
interproscan<- read.csv("~/iowa_state/lab/genome annotation files/ef2 blumeria genome/dh14_interproscan.csv", row.names = 1)
term2gene2=interproscan[, c("interproscan", "gene")]
term2name2=interproscan[, c("interproscan", "description")]
#uniqueGO<- read.csv("~/iowa_state/lab/GRN_reconstruction/HORVUR2 version analysis/ensemble of GRN/GO analyses/uniqueGO.csv")
unique_interproscan<- read.csv("~/iowa_state/lab/genome annotation files/ef2 blumeria genome/dh14_unique_interproscan.csv", row.names = 1)

unique(bgh_spread_tables$consensus)
type=c("bgh_additive_Mla6_Bln1", "bgh_symmetric", "bgh_masked", "bgh_suppression", "bgh_pseudo_masked", "bgh_positive", "bgh_negative",
       "bgh_unique_eff_Mla6", "bgh_unique_eff_Rar3", "bgh_equal_eff_Mla6_Rar3", "bgh_predominant_Mla6", "bgh_predominant_Rar3")
type %in% unique(bgh_spread_tables$consensus)

interproscan_result<- data.frame()

for (t in type){
  if(t %in% bgh_spread_tables$consensus){
    DE_gene_list<- bgh_spread_tables[bgh_spread_tables$consensus==t,"gene"]
    x_horvu <- enricher(DE_gene_list, TERM2GENE=term2gene2, TERM2NAME=term2name2)
    if(!is.null(x_horvu)){
      if(nrow(x_horvu)>0){
        x_horvu<- merge(as.data.frame(x_horvu), unique_interproscan, by.x="ID", by.y="interproscan")
        x_horvu$group<- t
        interproscan_result<- rbind(interproscan_result, x_horvu)
      }
    }
  }
}

interproscan_result<- interproscan_result[!duplicated(interproscan_result),]
write.csv(interproscan_result, "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/GO/R3/interproscan_epistasis_bgh.csv")

length(unique(interproscan_result$ID))
length(unique(interproscan_result$Description))


#go analysis per timepoint 

go_result_time<- data.frame()

for (time in c("t0","t16","t20","t24","t32","t48")){
  for (t in type){
    if(t %in% bgh_spread_tables[,time]){
    DE_gene_list<- bgh_spread_tables[bgh_spread_tables[,time]==t,"gene"]
    x_horvu <- enricher(DE_gene_list, TERM2GENE=term2gene, TERM2NAME=term2name)
    if(!is.null(x_horvu)){
      if(nrow(x_horvu)>0){
        x_horvu<- merge(as.data.frame(x_horvu), uniqueGO, by.x="ID", by.y="GO")
        x_horvu$type<- t
        x_horvu$group<- paste0(time, "_", t)
        x_horvu$time<- time
        go_result_time<- rbind(go_result_time, x_horvu)
      }}
  }}}

go_result_time<- go_result_time[!duplicated(go_result_time),]
write.csv(go_result_time, "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/GO/R3/GO_epistasis_bgh_timepoints.csv")

length(unique(go_result$ID))
length(unique(go_result$Description))


#plot barley
library(tidyverse)
library(RColorBrewer)
#colors<- brewer.pal(n = 6, name = "Dark2")
colors  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
names(colors)<- c("additive_Mla6_Bln1", "symmetric", "masked", "suppression", "pseudo_masked", "positive", "negative")

go_result_time<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/GO/R3/GO_epistasis_barley_timepoints.csv", row.names = 1)
go_result_time<- go_result_time[go_result_time$p.adjust<0.005 & 
                                  go_result_time$type %in% c("additive_Mla6_Bln1", "symmetric", "masked", "suppression", "pseudo_masked", "positive", "negative"),]
go_result_time_summary<- go_result_time %>% group_by(time, GO_term, type) %>% filter(p.adjust==min(p.adjust))
go_result_time_summary[go_result_time_summary$ID=="GO:0016709","GO_term"]<- "Cytochrome P450"
go_result_time_summary$Description_GO<- paste0(go_result_time_summary$GO_term, ":(",go_result_time_summary$ID, ")" )

go_result_time_summary$order<- ifelse(go_result_time_summary$type=="additive_Mla6_Bln1", 1, 
                                      ifelse(go_result_time_summary$type=="symmetric", 2,
                                             ifelse(go_result_time_summary$type=="masked", 3, 
                                                    ifelse(go_result_time_summary$type=="suppression", 4, 
                                                           ifelse(go_result_time_summary$type=="pseudo_masked", 5,
                                                                  ifelse(go_result_time_summary$type=="positive", 6, 7))))))

#go_result_time_summary$Description_GO = with(go_result_time_summary, reorder(Description_GO, -order))
go_result_time_summary$order2<- 10*(as.numeric(substring(go_result_time_summary$time, 2))+1)*go_result_time_summary$order
go_result_time_summary$Description_GO = with(go_result_time_summary, reorder(Description_GO, -order2))

pdf(paste0("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/GO/R3/GO_epistasis_time_barley.pdf"), width = 30, height = 40, fonts = "ArialMT", pointsize = 30)
ggplot(data=go_result_time_summary, aes(x=factor(time), y=Description_GO))  +  geom_raster(aes(fill=type)) + 
  theme(text = element_text(size=40), axis.text.x = element_text(angle = 90))+ 
  xlab("group") +ylab("Top GO Enriched Terms")+ scale_fill_manual(values = colors)#+ theme_grey(base_size = 40)

ggplot(data=go_result_time_summary[go_result_time_summary$p.adjust<0.0001,], aes(x=factor(time), y=Description_GO))  +  geom_raster(aes(fill=type)) + 
  theme(text = element_text(size=40), axis.text.x = element_text(angle = 90))+ 
  xlab("group") +ylab("Top GO Enriched Terms")+ scale_fill_manual(values = colors)#+ theme_grey(base_size = 40)

dev.off()


library(RColorBrewer)
#colors<- brewer.pal(n = 3, name = "Dark2")
colors  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
names(colors)<- c("additive_Mla6_Bln1", "symmetric", "masked", "suppression", "pseudo_masked", "positive", "negative")

go_result<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/GO/R3/GO_epistasis_barley.csv", row.names = 1)
go_result<- go_result[go_result$p.adjust<0.005 & 
                                  go_result$group %in% c("additive_Mla6_Bln1", "symmetric", "masked", "suppression", "pseudo_masked", "positive", "negative"),]
go_result_summary<- go_result %>% group_by(GO_term, group) %>% filter(p.adjust==min(p.adjust))

go_result_summary[go_result_summary$ID=="GO:0016709","GO_term"]<- "Cytochrome P450"
go_result_summary$Description_GO<- paste0(go_result_summary$GO_term, ":(",go_result_summary$ID, ")" )
go_result_summary$order<- ifelse(go_result_summary$group=="additive_Mla6_Bln1", 1, 
                                      ifelse(go_result_summary$group=="symmetric", 2,
                                             ifelse(go_result_summary$group=="masked", 3, 
                                                    ifelse(go_result_summary$group=="suppression", 4, 
                                                           ifelse(go_result_summary$group=="pseudo_masked", 5,
                                                                  ifelse(go_result_summary$group=="positive", 6, 7))))))

go_result_summary$Description_GO = with(go_result_summary, reorder(Description_GO, -order))


pdf(paste0("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/GO/R3/GO_epistasis_consensus_barley.pdf"), width = 30, height = 30, fonts = "ArialMT", pointsize = 30)
ggplot(data=go_result_summary, aes(x=factor(group), y=Description_GO)) + 
  theme(text = element_text(size=40)) +  geom_raster(aes(fill=group)) + 
  xlab("group") +ylab("Top GO Enriched Terms")+ scale_fill_manual(values = colors)#+ theme_grey(base_size = 40)
ggplot(data=go_result_time_summary[go_result_time_summary$p.adjust<0.0001,], aes(x=factor(time), y=Description_GO))  +  geom_raster(aes(fill=type)) + 
  theme(text = element_text(size=40), axis.text.x = element_text(angle = 90))+ 
  xlab("group") +ylab("Top GO Enriched Terms")+ scale_fill_manual(values = colors)#+ theme_grey(base_size = 40)

dev.off()




#bgh

library(tidyverse)
library(RColorBrewer)
#colors<- brewer.pal(n = 7, name = "Dark2")
colors  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
names(colors)<- c("bgh_additive_Mla6_Bln1", "bgh_symmetric", "bgh_masked", "bgh_suppression", "bgh_pseudo_masked", "bgh_positive", "bgh_negative")
go_result_time<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/GO/R3/GO_epistasis_bgh_timepoints.csv", row.names = 1)
go_result_time<- go_result_time[go_result_time$p.adjust<0.05 & 
                                  go_result_time$type %in% c("bgh_additive_Mla6_Bln1", "bgh_symmetric", "bgh_masked", "bgh_suppression", "bgh_pseudo_masked", "bgh_positive", "bgh_negative"),]
go_result_time_summary<- go_result_time %>% group_by(time, GO_term, type) %>% filter(p.adjust==min(p.adjust))
go_result_time_summary[go_result_time_summary$ID=="GO:0016705","GO_term"]<- "oxidoreductase activity"
go_result_time_summary$Description_GO<- paste0(go_result_time_summary$GO_term, ":(",go_result_time_summary$ID, ")" )

go_result_time_summary$order<- ifelse(go_result_time_summary$type=="bgh_additive_Mla6_Bln1", 1, 
                                      ifelse(go_result_time_summary$type=="bgh_symmetric", 2,
                                             ifelse(go_result_time_summary$type=="bgh_masked", 3, 
                                                    ifelse(go_result_time_summary$type=="bgh_suppression", 4, 
                                                           ifelse(go_result_time_summary$type=="bgh_pseudo_masked", 5,
                                                                  ifelse(go_result_time_summary$type=="bgh_positive", 6, 7))))))

#go_result_time_summary$Description_GO = with(go_result_time_summary, reorder(Description_GO, -order))
go_result_time_summary$order2<- 10*(as.numeric(substring(go_result_time_summary$time, 2))+1)*go_result_time_summary$order
go_result_time_summary$Description_GO = with(go_result_time_summary, reorder(Description_GO, -order2))

ggplot(data=go_result_time_summary, aes(x=factor(time), y=Description_GO))  +  geom_raster(aes(fill=type)) + 
  theme(text = element_text(size=10), axis.text.x = element_text(angle = 90))+ 
  xlab("group") +ylab("Top GO Enriched Terms")+ scale_fill_manual(values = colors)#+ theme_grey(base_size = 40)


pdf(paste0("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/GO/R3/GO_epistasis_time_bgh.pdf"), width = 30, height = 25, fonts = "ArialMT", pointsize = 30)
ggplot(data=go_result_time_summary, aes(x=factor(time), y=Description_GO))  +  geom_raster(aes(fill=type)) + 
  theme(text = element_text(size=40), axis.text.x = element_text(angle = 90))+ 
  xlab("group") +ylab("Top GO Enriched Terms")+ scale_fill_manual(values = colors)#+ theme_grey(base_size = 40)

ggplot(data=go_result_time_summary[go_result_time_summary$p.adjust<0.001,], aes(x=factor(time), y=Description_GO))  +  geom_raster(aes(fill=type)) + 
  theme(text = element_text(size=40), axis.text.x = element_text(angle = 90))+ 
  xlab("group") +ylab("Top GO Enriched Terms")+ scale_fill_manual(values = colors)#+ theme_grey(base_size = 40)

dev.off()


library(RColorBrewer)
#colors<- brewer.pal(n = 3, name = "Dark2")
colors  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
names(colors)<- c("bgh_additive_Mla6_Bln1", "bgh_symmetric", "bgh_masked", "bgh_suppression", "bgh_pseudo_masked", "bgh_positive", "bgh_negative")
go_result<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/GO/R3/GO_epistasis_bgh.csv", row.names = 1)
go_result<- go_result[go_result$p.adjust<0.05 & 
                        go_result$group %in% c("bgh_additive_Mla6_Bln1", "bgh_symmetric", "bgh_masked", "bgh_suppression", "bgh_pseudo_masked", "bgh_positive", "bgh_negative"),]
go_result_summary<- go_result %>% group_by(GO_term, group) %>% filter(p.adjust==min(p.adjust))

go_result_summary[go_result_summary$ID=="GO:0016620","GO_term"]<- "redox activity"
go_result_summary$Description_GO<- paste0(go_result_summary$GO_term, ":(",go_result_summary$ID, ")" )
go_result_summary$order<- ifelse(go_result_summary$group=="bgh_additive_Mla6_Bln1", 1, 
                                 ifelse(go_result_summary$group=="bgh_symmetric", 2,
                                        ifelse(go_result_summary$group=="bgh_masked", 3, 
                                               ifelse(go_result_summary$group=="bgh_suppression", 4, 
                                                      ifelse(go_result_summary$group=="bgh_pseudo_masked", 5,
                                                             ifelse(go_result_summary$group=="bgh_positive", 6, 7))))))

go_result_summary$Description_GO = with(go_result_summary, reorder(Description_GO, -order))


pdf(paste0("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/GO/R3/GO_epistasis_consensus_bgh.pdf"), width = 30, height = 20, fonts = "ArialMT", pointsize = 30)
ggplot(data=go_result_summary, aes(x=factor(group), y=Description_GO)) + 
  theme(text = element_text(size=40)) +  geom_raster(aes(fill=group)) + 
  xlab("group") +ylab("Top GO Enriched Terms")+ scale_fill_manual(values = colors)#+ theme_grey(base_size = 40)
ggplot(data=go_result_time_summary[go_result_time_summary$p.adjust<0.001,], aes(x=factor(time), y=Description_GO))  +  geom_raster(aes(fill=type)) + 
  theme(text = element_text(size=40), axis.text.x = element_text(angle = 90))+ 
  xlab("group") +ylab("Top GO Enriched Terms")+ scale_fill_manual(values = colors)#+ theme_grey(base_size = 40)

dev.off()



################################
#plot Mla6 rar3 
library(tidyverse)
library(RColorBrewer)
#colors<- brewer.pal(n = 6, name = "Dark2")
colors = c("#00AFBB", "#E7B800", "#FC4E07", "#52854C", "#7570B3")
names(colors)<- c("unique_eff_Mla6", "unique_eff_Rar3", "equal_eff_Mla6_Rar3", "predominant_Rar3", "predominant_Mla6")

go_result_time<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/GO/R3/GO_epistasis_barley_timepoints.csv", row.names = 1)
go_result_time<- go_result_time[go_result_time$p.adjust<0.005 & 
                                  go_result_time$type %in% c("unique_eff_Mla6", "unique_eff_Rar3", "equal_eff_Mla6_Rar3", "predominant_Rar3"),]
go_result_time_summary<- go_result_time %>% group_by(time, GO_term, type) %>% filter(p.adjust==min(p.adjust))
go_result_time_summary[go_result_time_summary$ID=="GO:0016709","GO_term"]<- "Cytochrome P450"
go_result_time_summary$Description_GO<- paste0(go_result_time_summary$GO_term, ":(",go_result_time_summary$ID, ")" )

go_result_time_summary$order<- ifelse(go_result_time_summary$type=="unique_eff_Mla6", 1, 
                                      ifelse(go_result_time_summary$type=="unique_eff_Rar3", 2,
                                             ifelse(go_result_time_summary$type=="equal_eff_Mla6_Rar3", 3, 
                                                    ifelse(go_result_time_summary$type=="predominant_Rar3", 4, 
                                                           ifelse(go_result_time_summary$type=="predominant_Mla6", 5,6)))))
#go_result_time_summary$Description_GO = with(go_result_time_summary, reorder(Description_GO, -order))
go_result_time_summary$order2<- 10*(as.numeric(substring(go_result_time_summary$time, 2))+1)*go_result_time_summary$order
go_result_time_summary$Description_GO = with(go_result_time_summary, reorder(Description_GO, -order2))

pdf(paste0("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/GO/R3/GO_mla6_rar3_time_barley.pdf"), width = 30, height = 40, fonts = "ArialMT", pointsize = 30)
ggplot(data=go_result_time_summary, aes(x=factor(time), y=Description_GO))  +  geom_raster(aes(fill=type)) + 
  theme(text = element_text(size=40), axis.text.x = element_text(angle = 90))+ 
  xlab("group") +ylab("Top GO Enriched Terms")+ scale_fill_manual(values = colors)#+ theme_grey(base_size = 40)

ggplot(data=go_result_time_summary[go_result_time_summary$p.adjust<0.001,], aes(x=factor(time), y=Description_GO))  +  geom_raster(aes(fill=type)) + 
  theme(text = element_text(size=40), axis.text.x = element_text(angle = 90))+ 
  xlab("group") +ylab("Top GO Enriched Terms")+ scale_fill_manual(values = colors)#+ theme_grey(base_size = 40)

dev.off()


library(RColorBrewer)
#colors<- brewer.pal(n = 3, name = "Dark2")
colors = c("#00AFBB", "#E7B800", "#FC4E07", "#52854C", "#7570B3")
names(colors)<- c("unique_eff_Mla6", "unique_eff_Rar3", "equal_eff_Mla6_Rar3", "predominant_Rar3", "predominant_Mla6")

go_result<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/GO/R3/GO_epistasis_barley.csv", row.names = 1)
go_result<- go_result[go_result$p.adjust<0.005 & 
                        go_result$group %in% c("unique_eff_Mla6", "unique_eff_Rar3", "equal_eff_Mla6_Rar3", "predominant_Rar3", "predominant_Mla6"),]
go_result_summary<- go_result %>% group_by(GO_term, group) %>% filter(p.adjust==min(p.adjust))

go_result_summary[go_result_summary$ID=="GO:0016709","GO_term"]<- "Cytochrome P450"
go_result_summary$Description_GO<- paste0(go_result_summary$GO_term, ":(",go_result_summary$ID, ")" )
go_result_summary$order<- ifelse(go_result_summary$group=="unique_eff_Mla6", 1, 
                                      ifelse(go_result_summary$group=="unique_eff_Rar3", 2,
                                             ifelse(go_result_summary$group=="equal_eff_Mla6_Rar3", 3, 
                                                    ifelse(go_result_summary$group=="predominant_Rar3", 4, 
                                                           ifelse(go_result_summary$group=="predominant_Mla6", 5,6)))))

go_result_summary$Description_GO = with(go_result_summary, reorder(Description_GO, -order))


pdf(paste0("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/GO/R3/GO_mla6_rar3_consensus_barley.pdf"), width = 30, height = 30, fonts = "ArialMT", pointsize = 30)
ggplot(data=go_result_summary, aes(x=factor(group), y=Description_GO)) + 
  theme(text = element_text(size=40)) +  geom_raster(aes(fill=group)) + 
  xlab("group") +ylab("Top GO Enriched Terms")+ scale_fill_manual(values = colors)#+ theme_grey(base_size = 40)
ggplot(data=go_result_time_summary[go_result_time_summary$p.adjust<0.001,], aes(x=factor(time), y=Description_GO))  +  geom_raster(aes(fill=type)) + 
  theme(text = element_text(size=40), axis.text.x = element_text(angle = 90))+ 
  xlab("group") +ylab("Top GO Enriched Terms")+ scale_fill_manual(values = colors)#+ theme_grey(base_size = 40)

dev.off()




#bgh

library(tidyverse)
library(RColorBrewer)
#colors<- brewer.pal(n = 7, name = "Dark2")
colors = c("#00AFBB", "#E7B800", "#FC4E07", "#52854C", "#7570B3")
names(colors)<- c("bgh_unique_eff_Mla6", "bgh_unique_eff_Rar3", "bgh_equal_eff_Mla6_Rar3", "bgh_predominant_Mla6", "bgh_predominant_Rar3")

go_result_time<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/GO/R3/GO_epistasis_bgh_timepoints.csv", row.names = 1)
go_result_time<- go_result_time[go_result_time$p.adjust<0.05 & 
                                  go_result_time$type %in% c("bgh_unique_eff_Mla6", "bgh_unique_eff_Rar3", "bgh_equal_eff_Mla6_Rar3", "bgh_predominant_Rar3", "bgh_predominant_Mla6"),]
go_result_time_summary<- go_result_time %>% group_by(time, GO_term, type) %>% filter(p.adjust==min(p.adjust))
go_result_time_summary[go_result_time_summary$ID=="GO:0016705","GO_term"]<- "oxidoreductase activity"
go_result_time_summary$Description_GO<- paste0(go_result_time_summary$GO_term, ":(",go_result_time_summary$ID, ")" )

go_result_time_summary$order<- ifelse(go_result_time_summary$type=="bgh_unique_eff_Mla6", 1, 
                                      ifelse(go_result_time_summary$type=="bgh_unique_eff_Rar3", 2,
                                             ifelse(go_result_time_summary$type=="bgh_equal_eff_Mla6_Rar3", 3, 
                                                    ifelse(go_result_time_summary$type=="bgh_predominant_Rar3", 4, 
                                                           ifelse(go_result_time_summary$type=="bgh_predominant_Mla6", 5, 6)))))

#go_result_time_summary$Description_GO = with(go_result_time_summary, reorder(Description_GO, -order))
go_result_time_summary$order2<- 10*(as.numeric(substring(go_result_time_summary$time, 2))+1)*go_result_time_summary$order
go_result_time_summary$Description_GO = with(go_result_time_summary, reorder(Description_GO, -order2))

ggplot(data=go_result_time_summary, aes(x=factor(time), y=Description_GO))  +  geom_raster(aes(fill=type)) + 
  theme(text = element_text(size=10), axis.text.x = element_text(angle = 90))+ 
  xlab("group") +ylab("Top GO Enriched Terms")+ scale_fill_manual(values = colors)#+ theme_grey(base_size = 40)


pdf(paste0("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/GO/R3/GO_mla6_rar3_time_bgh.pdf"), width = 30, height = 25, fonts = "ArialMT", pointsize = 30)
ggplot(data=go_result_time_summary, aes(x=factor(time), y=Description_GO))  +  geom_raster(aes(fill=type)) + 
  theme(text = element_text(size=40), axis.text.x = element_text(angle = 90))+ 
  xlab("group") +ylab("Top GO Enriched Terms")+ scale_fill_manual(values = colors)#+ theme_grey(base_size = 40)

ggplot(data=go_result_time_summary[go_result_time_summary$p.adjust<0.001,], aes(x=factor(time), y=Description_GO))  +  geom_raster(aes(fill=type)) + 
  theme(text = element_text(size=40), axis.text.x = element_text(angle = 90))+ 
  xlab("group") +ylab("Top GO Enriched Terms")+ scale_fill_manual(values = colors)#+ theme_grey(base_size = 40)

dev.off()


library(RColorBrewer)
#colors<- brewer.pal(n = 3, name = "Dark2")
colors = c("#00AFBB", "#E7B800", "#FC4E07", "#52854C", "#7570B3")
names(colors)<- c("bgh_unique_eff_Mla6", "bgh_unique_eff_Rar3", "bgh_equal_eff_Mla6_Rar3", "bgh_predominant_Mla6", "bgh_predominant_Rar3")
go_result<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/GO/R3/GO_epistasis_bgh.csv", row.names = 1)
go_result<- go_result[go_result$p.adjust<0.05 & 
                        go_result$group %in% c("bgh_unique_eff_Mla6", "bgh_unique_eff_Rar3", "bgh_equal_eff_Mla6_Rar3", "bgh_predominant_Rar3", "bgh_predominant_Mla6"),]
go_result_summary<- go_result %>% group_by(GO_term, group) %>% filter(p.adjust==min(p.adjust))

go_result_summary[go_result_summary$ID=="GO:0016620","GO_term"]<- "redox activity"
go_result_summary$Description_GO<- paste0(go_result_summary$GO_term, ":(",go_result_summary$ID, ")" )
go_result_summary$order<- ifelse(go_result_summary$group=="bgh_unique_eff_Mla6", 1, 
                                 ifelse(go_result_summary$group=="bgh_unique_eff_Rar3", 2,
                                        ifelse(go_result_summary$group=="bgh_equal_eff_Mla6_Rar3", 3, 
                                               ifelse(go_result_summary$group=="bgh_predominant_Rar3", 4, 
                                                      ifelse(go_result_summary$group=="bgh_predominant_Mla6", 5, 6)))))

go_result_summary$Description_GO = with(go_result_summary, reorder(Description_GO, -order))


pdf(paste0("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/GO/R3/GO_mla6_rar3_consensus_bgh.pdf"), width = 30, height = 20, fonts = "ArialMT", pointsize = 30)
ggplot(data=go_result_summary, aes(x=factor(group), y=Description_GO)) + 
  theme(text = element_text(size=40)) +  geom_raster(aes(fill=group)) + 
  xlab("group") +ylab("Top GO Enriched Terms")+ scale_fill_manual(values = colors)#+ theme_grey(base_size = 40)
ggplot(data=go_result_time_summary[go_result_time_summary$p.adjust<0.001,], aes(x=factor(time), y=Description_GO))  +  geom_raster(aes(fill=type)) + 
  theme(text = element_text(size=40), axis.text.x = element_text(angle = 90))+ 
  xlab("group") +ylab("Top GO Enriched Terms")+ scale_fill_manual(values = colors)#+ theme_grey(base_size = 40)

dev.off()






