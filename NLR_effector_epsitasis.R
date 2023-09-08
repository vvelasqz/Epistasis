# characterize NLRs and effectors under epistasis models

expression <- read.csv("~/iowa_state/lab/RNAseq/DESeq2 HvR3/hv_R3_genes_norm_counts_de_tax_sp.csv", stringsAsFactors = F)
spread_tables_Mla6_Bln1<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/spread_tables_Mla6_Bln1.csv", row.names = 1)
spread_tables_Mla6_Rar3<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/spread_tables_Mla6_Rar3.csv", row.names = 1)
spread_tables_barley<- rbind(spread_tables_Mla6_Bln1, spread_tables_Mla6_Rar3)

NLR<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/R3_NLRs.csv")
sum(NLR$gene %in% spread_tables_Mla6_Bln1$gene)
sum(NLR$gene %in% spread_tables_Mla6_Rar3$gene)
sum(NLR$gene %in% expression$gene)

mla_interactors_r3<- read.csv("~/iowa_state/lab/Y2H-Seq effectors/results/process_output_Y2H-SCORES/mla_interactors_r3.csv", row.names = 1)
sum(NLR$gene %in% mla_interactors_r3)
interactome_R2_HC<- read.csv(file='~/iowa_state/lab/interactome/combine data spp/tritex interactome V2/d_weighted_interactome_R2_HC2.csv')
convert_barley_V2_R2_R3<- read.csv("~/iowa_state/lab/genome annotation files/tritex R3/convert_barley_V2_R2_R3.csv")[,c(1,10)]
interactome_R3_HC<- merge(interactome_R2_HC, convert_barley_V2_R2_R3, by.x="source", by.y = "gene_R2", all.x = T)
interactome_R3_HC<- merge(interactome_R3_HC, convert_barley_V2_R2_R3, by.x="target", by.y = "gene_R2", all.x=T)
interactome_R3_HC<- interactome_R3_HC[!duplicated(interactome_R3_HC),][,c(4,5,3)]
colnames(interactome_R3_HC)[1:2]<- c("source", "target")
nodes_hvint<- unique(c(interactome_R3_HC$source,interactome_R3_HC$target ))
sum(NLR$gene %in% nodes_hvint)
NLR_interactome<- interactome_R3_HC[interactome_R3_HC$source %in% NLR$gene | interactome_R3_HC$target %in% NLR$gene,]

epi_NLR<- merge(NLR[,2:3], spread_tables_barley, by="gene")
epi_NLR_interactome<- interactome_R3_HC[interactome_R3_HC$source %in% epi_NLR$gene | interactome_R3_HC$target %in% epi_NLR$gene,]
epi_NLR_interactome <- merge(interactome_R3_HC, spread_tables_barley[,c(1,8)], by.x="source", by.y ="gene", all.x = T)
epi_NLR_interactome <- merge(epi_NLR_interactome, spread_tables_barley[,c(1,8)], by.x="target", by.y ="gene", all.x = T)
epi_NLR_interactome<- epi_NLR_interactome[!duplicated(epi_NLR_interactome),]
epi_NLR_interactome<- epi_NLR_interactome[epi_NLR_interactome$source %in% epi_NLR$gene | epi_NLR_interactome$target %in% epi_NLR$gene,]
sum(epi_NLR_interactome$source %in% epi_NLR$gene)
sum(epi_NLR_interactome$target %in% epi_NLR$gene)


effectors<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/effectors.csv")
bgh_spread_tables_Mla6_Bln1<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/bgh/HVR3/bgh_spread_tables_Mla6_Bln1.csv", row.names = 1)
bgh_spread_tables_Mla6_Rar3<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/bgh/HVR3/bgh_spread_tables_Mla6_Rar3.csv", row.names = 1)
sum(effectors$EF2_gene_ID %in% bgh_spread_tables_Mla6_Bln1$gene)
sum(effectors$EF2_gene_ID %in%  bgh_spread_tables_Mla6_Rar3$gene)
sum(effectors$EF2_gene_ID %in% expression$gene)

spread_tables_bgh<- rbind(bgh_spread_tables_Mla6_Bln1, bgh_spread_tables_Mla6_Rar3)
epi_effectors<- merge(spread_tables_bgh, effectors, by.y="EF2_gene_ID",by.x="gene")

write.csv(epi_NLR, "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/epi_NLR.csv")
write.csv(epi_effectors, "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/epi_effectors.csv")

#plot NLR and effector expression
library(gplots)
library(multcompView)
library(reshape2)
#taxon specific normalization
#df = read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/small_RNAs/Bgh_sRNA_read_counts_filtered.csv", row.names = 1)
df = read.csv("~/iowa_state/lab/RNAseq/DESeq2 HvR3/hv_R3_genes_norm_counts_de_tax_sp.csv")
hvu_de = read.csv("~/iowa_state/lab/RNAseq/DESeq2 HvR3/hv_R3_genes_deseq2_results_pairwise_genotype_padj_tax_sp.csv")

#barley#annotation barley or blumeria
#Significance threshold for p-adj, barley 0.001 or bhg 0.05
thr <- 0.001
ann = read.csv("~/iowa_state/lab/genome annotation files/tritex R3/HR3_annotation.csv",row.names = 1)[,c(2,5)]
colnames(ann)<- c("gene", "description")
#if plotting all NLRs
preys <- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/R3_NLRs.csv")
preys$group<- "NLR"
chrs = unique(preys$group)

#if plotting subgroups consensus
preys <- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/epi_NLR.csv", row.names = 1)
preys<- preys[preys$consensus!="",]
preys$group<- preys$consensus
chrs = unique(preys$consensus)


#if plotting subgroups time
preys <- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/epi_NLR.csv", row.names = 1)
preys<- melt(preys,measure.vars = 3:8)
preys<- preys[!is.na(preys$value),]
#preys[is.na(preys)]<-"NA"
#chrs = unique(preys$consensus)
preys$group<- paste0(preys$variable,"_",preys$value)
chrs = unique(preys$group)

#blumeria
thr <- 0.003
ann = read.table("~/iowa_state/lab/RNAseq/DESeq2 analysis HvV2/bgh_ef1_ef2_CSEP.tsv", header = T)[,c(2,11)]
colnames(ann)<- c("gene","description")
preys <- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/epi_effectors.csv", row.names = 1)
preys$group<- preys$consensus
#preys[is.na(preys)]<-"NA"
chrs = unique(preys$consensus)
chrs<- chrs[!chrs %in% ""]

#if plotting all epi effectors or all effectors
#preys <- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/epi_effectors.csv", row.names = 1)

preys<- data.frame(gene=read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/effectors.csv")[,9])

length(unique(preys$gene))
preys$group<- "effector"
chrs = unique(preys$group)

  
# df = read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/small_RNAs/Bgh_sRNA_read_counts_normalized.csv")
# colnames(df)[1]<- "gene"
# hvu_de = read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/small_RNAs/srna_genes_deseq2_results_pairwise_genotype_padj.csv")
# hvu_de[is.na(hvu_de)]<- 1

#Differential expression results
hvu_de[is.na(hvu_de)]<-1
#Put it all together
df <- subset(df, df$gene %in% grep(paste(unique(preys$gene), collapse="|"), df$gene, value = TRUE))
data = merge(ann, df, all.y = T)
data = merge(data, hvu_de, all.x = TRUE)
data[is.na(data)]<- 1

for(m in 1:length(chrs)){
  dc = subset(data, data$gene %in% preys$gene[preys$group==chrs[m]])
  #dc = subset(data, data$gene %in% preys$gene[preys$consensus==chrs[m]])  # select the preys that correspond to the bait m
  #dc = subset(data, data$gene %in% grep(chrs[m], data$gene, value = TRUE))
  n = nrow(dc)
  print(n)
  
  #filename = paste0("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/graphs_consensus/graphs_NLRs_sig_letters_", chrs[m], ".pdf")
  filename = paste0("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/graphs_All_effector_sig_letters_", chrs[m], ".pdf")
  #filename = paste0("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/graphs_effector_sig_letters_", chrs[m], ".pdf")
  #filename = paste0("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/graphs_consensus/graphs_effector_sig_letters_", chrs[m], ".pdf")
  #filename = paste0("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/graphs_NLRs_sig_letters_", chrs[m], ".pdf")
  #filename = paste0("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/graphs_time/graphs_NLRs_sig_letters_", chrs[m], ".pdf")
  #filename = paste0("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/graphs_time/graphs_effector_sig_letters_", chrs[m], ".pdf")
  pdf(filename, width = 10, height = 10, fonts = "ArialMT", pointsize = 12)
  
  for(i in 1:n){
    # I will calculate the significance groups across genotype/timepoint comparisons
    symbs2 = dc[i, 93:152]
    letterCodes <- data.frame()
    for (time in paste0("t", c(0, 16, 20, 24, 32, 48))) {
      symbs2_time<- symbs2[,c(grep(time, colnames(symbs2), value = TRUE))]
      vector_multC <- as.numeric(symbs2_time[1,])
      names(vector_multC) <- apply(data.frame(colnames(symbs2_time)), 1, FUN= function(x){
        paste0(strsplit(x, "_")[[1]][2], "-", strsplit(x, "_")[[1]][3])})
      compLetters <- multcompLetters(vector_multC, threshold = thr)
      if(ncol(letterCodes)==0){
        letterCodes <- data.frame(compLetters$Letters, stringsAsFactors=FALSE)
      }else{
        letterCodes <- cbind(letterCodes, data.frame(compLetters$Letters, stringsAsFactors=FALSE))}
      colnames(letterCodes)[ncol(letterCodes)]<- time
      rm(vector_multC, compLetters)
    }
    
    gi = dc[i,3:92]
    
    # Calculating means for gi
    gmeans = c()
    j = 1
    for(k in 1:30){
      gmeans[k] = mean(as.numeric(gi[1,j:(j+2)]))
      j = j+3
    }
    rm(j, k, gi)
    
    # Reading in means for genotype, depending on the order in the table
    
    # For V2 new annotation 
    means.bln1.norm = gmeans[19:24]
    means.dm.norm = gmeans[25:30]
    means.mla6.norm = gmeans[7:12]
    means.rar3.norm = gmeans[13:18]
    means.wt.norm = gmeans[1:6]
    
    # Calculating y limit
    down2=min(gmeans[1:30])
    if(is.na(down2) == TRUE){down2 = 0}
    up2=max(gmeans[1:30])
    if(is.na(up2) == TRUE){up2 = down2}
    lim2=c(down2, up2*1.1)
    
    # Creating a nice title
    title = paste(sapply(dc[i,c("gene","description")], as.character), collapse = "  |  ") 
    
    # Plotting
    plotCI(c(1,5,6,7,9,13), unlist(means.mla6.norm), lty=2, xaxt="n", xlim=c(0.5,13.5), ylim=as.numeric(lim2), gap=0, ylab="DESeq2 Normalized Read Counts by Timepoint", xlab="Hours after Inoculation", main=title, type="n", col="red")
    lines(c(1,5,6,7,9,13), unlist(means.mla6.norm), col="red", lty=2)
    text(c(1,5,6,7,9,13), unlist(means.mla6.norm), as.character(letterCodes["mla6",]), col="red")
    
    plotCI(c(0.8,4.8,5.8,6.8,8.8,12.8), unlist(means.wt.norm), lty=1, xaxt="n", xlim=c(0.3,13.7), ylim=as.numeric(lim2), gap=0, add=TRUE, type="n", col="black")
    lines(c(0.8,4.8,5.8,6.8,8.8,12.8), unlist(means.wt.norm), col="black", lty=1)
    text(c(0.8,4.8,5.8,6.8,8.8,12.8), unlist(means.wt.norm), as.character(letterCodes["wt",]), col="black")
    
    plotCI(c(0.9,4.9,5.9,6.9,8.9,12.9), unlist(means.bln1.norm), lty=1, xaxt="n", xlim=c(0.3,13.7), ylim=as.numeric(lim2), gap=0, add=TRUE, type="n", col="brown")
    lines(c(0.9,4.9,5.9,6.9,8.9,12.9), unlist(means.bln1.norm), col="brown", lty=1)
    text(c(0.9,4.9,5.9,6.9,8.9,12.9),  unlist(means.bln1.norm), as.character(letterCodes["bln1",]), col="brown")
    
    plotCI(c(1.1,5.1,6.1,7.1,9.1,13.1), unlist(means.rar3.norm), lty=2, xaxt="n", xlim=c(0.3,13.7), ylim=as.numeric(lim2), gap=0, add=TRUE, type="n", col="darkgreen")
    lines(c(1.1,5.1,6.1,7.1,9.1,13.1), unlist(means.rar3.norm), col="darkgreen", lty=2)
    text(c(1.1,5.1,6.1,7.1,9.1,13.1), unlist(means.rar3.norm), as.character(letterCodes["rar3",]), col="darkgreen")
    
    plotCI(c(1.2,5.2,6.2,7.2,9.2,13.2), unlist(means.dm.norm), lty=2, xaxt="n", xlim=c(0.3,13.7), ylim=as.numeric(lim2), gap=0, add=TRUE, type="n", col="blue")
    lines(c(1.2,5.2,6.2,7.2,9.2,13.2), unlist(means.dm.norm), col="blue", lty=2)
    text(c(1.2,5.2,6.2,7.2,9.2,13.2), unlist(means.dm.norm), as.character(letterCodes["dm",]), col="blue")
    
    axis(side=1, at=c(1,5,6,7,9,13), labels=c(0,16,20,24,32,48), cex=0.7)
    legend("top", legend=c("CI 16151", expression(italic("bln1-m19089")), expression(italic("mla6-m18982")), expression(italic("rar3-m11526")), expression(italic("(mla6+bln1)-m19028"))),lty=c(1,1,2,2,2), col=c("black","brown","red","darkgreen","blue"), bty = "n", horiz = TRUE, text.col = c("black","brown","red","darkgreen","blue"), cex = 0.75)
    
    rm(gmeans, lim2)
    
  }
  
  dev.off()
}


#MLA HORVU.MOREX.r3.1HG0012670





#plot whole matrix
library(tidyverse)
library(plotrix)
library(viridis)
source("~/iowa_state/lab/MLA GRN publication/tritex V2/interactome_analysis/color_matrix.R")
epi_NLR<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/epi_NLR.csv", row.names = 1)
rw_table<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/rw_types.csv")
epi_NLR<- merge(epi_NLR, rw_table, by="gene",all.x = T)

epi_effectors<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/epi_effectors.csv", row.names = 1)

#reorder
# counts<- counts[, c(grep("wt", colnames(counts)),
#                     grep("mla6", colnames(counts)),
#                     grep("rar3", colnames(counts)),
#                     grep("bln1", colnames(counts)),
#                     grep("dm", colnames(counts)))]
# colnames(counts)
#in the order of the data in counts
#gen=unique(sapply(colnames(counts), FUN=function(x) strsplit(x, "_")[[1]][1]))
#gen=c("wt", "bln1", "dm","mla6", "rar3")
#colnames(counts) <- paste0(rep(gen, each=3*6), "_t",rep(rep(time, each=3),5), "_r",rep(c(1,2,3),5*6))
time=c(0, 16, 20, 24, 32, 48)
gen=c("wt", "mla6", "rar3","bln1", "dm")
expression <- read.csv("~/iowa_state/lab/RNAseq/DESeq2 HvR3/hv_R3_genes_norm_counts_de_tax_sp.csv", stringsAsFactors = F, row.names = 1)
expression_mean<- matrix(nrow = nrow(expression), ncol = 30, dimnames = list(rownames(expression), paste0(rep(gen, each=6), "_t",rep(time,5))))
for(i in 1:nrow(expression)){
  j = 1
  for(k in 1:30){
    expression_mean[i,k] = mean(as.numeric(expression[i,j:(j+2)]))
    j = j+3
  }
}

dim(expression_mean)

#reorder
expression_mean<- expression_mean[, c(grep("rar3", colnames(expression_mean)),
                    grep("wt", colnames(expression_mean)),
                    grep("mla6", colnames(expression_mean)),
                    grep("dm", colnames(expression_mean)),
                    grep("bln1", colnames(expression_mean)))]
colnames(expression_mean)

write.csv(expression_mean, "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/expression_mean.csv")
expression_mean<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/expression_mean.csv", row.names = 1)

patterns<- c("additive_Mla6_Bln1", "symmetric", "masked", "suppression", "pseudo_masked", "positive", "negative", NA)
epi_NLR_mla6_bln1<- epi_NLR[epi_NLR$consensus %in% patterns, ]
colors<- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#FFFFFF")
names(colors)<- patterns
epi_NLR_mla6_bln1$color<- sapply(epi_NLR_mla6_bln1$consensus, FUN=function(x)colors[x] )

NLR_exp<- log(as.matrix(expression_mean[grep(paste(unique(epi_NLR$gene), collapse="|"), rownames(expression_mean)), ])+1)
length(unique(epi_NLR$gene))
#NLR_exp<- as.matrix(expression_mean[grep(paste(unique(epi_NLR$gene), collapse="|"), rownames(expression_mean)), ])

NLR_mla6_bln1_exp<- log(as.matrix(expression_mean[grep(paste(unique(epi_NLR_mla6_bln1$gene), collapse="|"), rownames(expression_mean)), ])+1)
ordered_heatmap<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/ordered_heatmap.csv")
epi_NLR_mla6_bln1<- epi_NLR_mla6_bln1[match(ordered_heatmap$ordered_heatmap, epi_NLR_mla6_bln1$gene),]
NLR_mla6_bln1_exp<- NLR_mla6_bln1_exp[match(epi_NLR_mla6_bln1$gene, row.names(NLR_mla6_bln1_exp)),]
table(epi_NLR_mla6_bln1$consensus)
col_heatmap<- unique(c(NLR_exp))
names(col_heatmap)<-viridis(length(unique(c(NLR_exp))))


pdf(paste0("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/heatmaps_mla6_bln1.pdf"), width = 50, height = 30, fonts = "ArialMT", pointsize = 30)
heatmap(NLR_mla6_bln1_exp, Colv=NA, col=names(col_heatmap[col_heatmap %in% unique(NLR_exp[rownames(NLR_exp) %in% epi_NLR_mla6_bln1$gene])]), RowSideColors=epi_NLR_mla6_bln1$color)
for(pattern in c("additive_Mla6_Bln1", "symmetric", "masked", "suppression", "negative")){
  #pattern<-"masked"
  heatmap(NLR_mla6_bln1_exp[rownames(NLR_mla6_bln1_exp) %in% epi_NLR_mla6_bln1[epi_NLR_mla6_bln1$consensus==pattern,"gene"],],
          Colv=NA, col= names(col_heatmap[col_heatmap %in% unique(NLR_exp[rownames(NLR_exp) %in% epi_NLR_mla6_bln1[epi_NLR_mla6_bln1$consensus==pattern,"gene"]])]), RowSideColors=epi_NLR_mla6_bln1[epi_NLR_mla6_bln1$consensus==pattern,"color"])
}
heatmap(NLR_mla6_bln1_exp, Rowv = NA, Colv=NA, col=names(col_heatmap[col_heatmap %in% unique(NLR_exp[rownames(NLR_exp) %in% epi_NLR_mla6_bln1$gene])]), RowSideColors=epi_NLR_mla6_bln1$color)
#heatmap(NLR_mla6_bln1_exp, Rowv = NA, Colv=NA, col=names(col_heatmap[col_heatmap %in% unique(NLR_exp[rownames(NLR_exp) %in% epi_NLR_mla6_bln1$gene])]), RowSideColors=epi_NLR_mla6_bln1$color)
dev.off()

#heatmap with everything annotating for mla6 rar3 and mla6 bln1
library(tidyverse)
require(scales)
library(munsell)
library(viridis)
library(gprofiler2)
#wt_log_scaled$`0` <- seq_gradient_pal("darkblue", "yellow")(wt_log_scaled$`0`)

epi_NLR<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/epi_NLR.csv", row.names = 1)
expression_mean<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/expression_mean.csv", row.names = 1)
patterns<- c("additive_Mla6_Bln1", "symmetric", "masked", "suppression", "pseudo_masked", "positive", "negative", NA)
epi_NLR_mla6_bln1<- epi_NLR[epi_NLR$consensus %in% patterns, ]
colors<- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#FFFFFF")
names(colors)<- patterns
epi_NLR_mla6_bln1$color<- sapply(epi_NLR_mla6_bln1$consensus, FUN=function(x)colors[x] )
patterns_rar3<- c( "unique_eff_Mla6", "unique_eff_Rar3", "equal_eff_Mla6_Rar3", "predominant_Mla6", "predominant_Rar3", NA)
epi_NLR_mla6_rar3<- epi_NLR[epi_NLR$consensus %in% patterns_rar3, ]
epi_NLR2<- merge(epi_NLR_mla6_bln1[,c(1,9)], epi_NLR_mla6_rar3[,c(1,9)], by="gene", all=T)
ordered_heatmap2<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/ordered_heatmap_complete.csv")

colors_rar3<-  c("#00AFBB", "#E7B800", "#FC4E07", "#52854C", "#7570B3", "#FFFFFF")
names(colors_rar3)<- patterns_rar3
epi_NLR2$color_bln1<- sapply(epi_NLR2$consensus.x, FUN=function(x)colors[x] )
epi_NLR2$color_rar3<- sapply(epi_NLR2$consensus.y, FUN=function(x)colors_rar3[x] )
epi_NLR2<- epi_NLR2[match(ordered_heatmap2$ordered_heatmap, epi_NLR2$gene),]

NLR_exp<- log(as.matrix(expression_mean[grep(paste(unique(epi_NLR$gene), collapse="|"), rownames(expression_mean)), ])+1)
# length(unique(epi_NLR$gene))
# epi_NLR2<- spread(epi_NLR[,c(1,9)]) #, key = "type", value = "type")

NLR_exp_filtered<- NLR_exp[rownames(NLR_exp) %in% epi_NLR2$gene,]
NLR_exp_filtered<- NLR_exp_filtered[match(ordered_heatmap2$ordered_heatmap, row.names(NLR_exp_filtered)),]

pdf(paste0("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/heatmaps_all.pdf"), width = 50, height = 30, fonts = "ArialMT", pointsize = 30)
heatmap(NLR_exp_filtered, Rowv = NA, Colv=NA, col=viridis(length(unique(c(NLR_exp_filtered)))), RowSideColors=epi_NLR2$color_bln1)
heatmap(NLR_exp_filtered, Rowv = NA, Colv=NA, col=viridis(length(unique(c(NLR_exp_filtered)))), RowSideColors=epi_NLR2$color_rar3)
dev.off()

cols<- sort(unique(c(NLR_exp_filtered)))
hist(cols)
cols_scaled <- ((cols-min(cols, na.rm=TRUE))/(max(cols, na.rm=TRUE)- min(cols, na.rm=TRUE)))
hist(cols_scaled)

mapViridis(cols, domain_min = min(cols), domain_max = max(cols), n = 256)
length(mapViridis(cols, domain_min = min(cols), domain_max = max(cols), n = 256))
length(viridis(length(unique(c(NLR_exp_filtered)))))

pdf(paste0("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/heatmaps_all2.pdf"), width = 50, height = 30, fonts = "ArialMT", pointsize = 30)
heatmap(NLR_exp_filtered, Rowv = NA, Colv=NA, col=viridis(length(unique(c(NLR_exp_filtered)))), RowSideColors=epi_NLR2$color_bln1)
heatmap(NLR_exp_filtered, Rowv = NA, Colv=NA, col=viridis(length(unique(c(NLR_exp_filtered))), option = "F"), RowSideColors=epi_NLR2$color_bln1)
heatmap(NLR_exp_filtered, Rowv = NA, Colv=NA, col=seq_gradient_pal("darkblue", "yellow")(cols_scaled), RowSideColors=epi_NLR2$color_rar3)
#heatmap(NLR_exp_filtered, Rowv = NA, Colv=NA, col=mapViridis(cols, domain_min = min(cols), domain_max = max(cols), n = length(cols)), RowSideColors=epi_NLR2$color_rar3)
dev.off()

pdf(paste0("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/heatmaps_all_clustered.pdf"), width = 50, height = 30, fonts = "ArialMT", pointsize = 30)
heatmap(NLR_exp_filtered, Colv=NA, col=viridis(length(unique(c(NLR_exp_filtered)))), RowSideColors=epi_NLR2$color_bln1)
heatmap(NLR_exp_filtered, Colv=NA, col=viridis(length(unique(c(NLR_exp_filtered)))), RowSideColors=epi_NLR2$color_rar3)
dev.off()


#now effector heatmap
epi_eff<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/epi_effectors.csv", row.names = 1)[,c(1,8)]
expression_mean<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/expression_mean.csv", row.names = 1)
unique(epi_effectors$consensus)

patterns<- c("bgh_additive_Mla6_Bln1", "bgh_symmetric", "bgh_masked", "bgh_suppression", "bgh_pseudo_masked", "bgh_positive", "bgh_negative", NA, "")
epi_eff_mla6_bln1<- epi_eff[epi_eff$consensus %in% patterns[1:7], ]
colors<- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#FFFFFF", "#FFFFFF")
names(colors)<- patterns
patterns_rar3<- c( "bgh_unique_eff_Mla6", "bgh_unique_eff_Rar3", "bgh_equal_eff_Mla6_Rar3", "bgh_predominant_Mla6", "bgh_predominant_Rar3", NA, "")
epi_eff_mla6_rar3<- epi_eff[epi_eff$consensus %in% patterns_rar3[1:5], ]
epi_eff2<- merge(epi_eff_mla6_bln1, epi_eff_mla6_rar3, by="gene", all=T)

colors_rar3<-  c("#00AFBB", "#E7B800", "#FC4E07", "#52854C", "#7570B3", "#FFFFFF", "#FFFFFF")
names(colors_rar3)<- patterns_rar3
epi_eff2$color_bln1<- sapply(epi_eff2$consensus.x, FUN=function(x)colors[x] )
epi_eff2$color_rar3<- sapply(epi_eff2$consensus.y, FUN=function(x)colors_rar3[x] )

eff_exp<- log(as.matrix(expression_mean[grep(paste(unique(epi_eff2$gene), collapse="|"), rownames(expression_mean)), ])+1)
length(unique(epi_eff$gene))

epi_eff2<- epi_eff2[!duplicated(epi_eff2$gene) & epi_eff2$gene %in% rownames(eff_exp),]


eff_mla6_bln1_exp<- log(as.matrix(expression_mean[grep(paste(unique(epi_eff_mla6_bln1$gene), collapse="|"), rownames(expression_mean)), ])+1)

heatmap(eff_mla6_bln1_exp, Colv=NA, col=viridis(length(unique(c(eff_exp)))), RowSideColors=epi_eff2[epi_eff2$gene %in% rownames(eff_mla6_bln1_exp),"color_bln1"])

pdf(paste0("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/heatmaps_effectors_clustered_long.pdf"), width = 50, height = 70, fonts = "ArialMT", pointsize = 30)
heatmap(eff_exp, Colv=NA, col=viridis(length(unique(c(eff_exp)))), RowSideColors=epi_eff2[epi_eff2$gene %in% rownames(eff_exp),"color_bln1"])
heatmap(eff_exp, Colv=NA, col=viridis(length(unique(c(eff_exp)))), RowSideColors=epi_eff2[epi_eff2$gene %in% rownames(eff_exp),"color_rar3"])
dev.off()


#make mini graphs

# Modified by Valeria to select preys in relationship with the respective bait interactors highlighting the DE points mini
#loading libraries
library(gplots)
library(multcompView)
#set working directory
df = read.csv("~/iowa_state/lab/RNAseq/DESeq2 HvR3/hv_R3_genes_norm_counts_de_tax_sp.csv")
hvu_de = read.csv("~/iowa_state/lab/RNAseq/DESeq2 HvR3/hv_R3_genes_deseq2_results_pairwise_genotype_padj_tax_sp.csv")

#barley#annotation barley or blumeria
#Significance threshold for p-adj, barley 0.001 or bhg 0.05
thr <- 0.001
ann = read.csv("~/iowa_state/lab/genome annotation files/tritex R3/HR3_annotation.csv",row.names = 1)[,c(2,5)]
colnames(ann)<- c("gene", "description")
#if plotting all NLRs
preys <- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/R3_NLRs.csv")
preys$group<- "NLR"
chrs = unique(preys$group)

#if plotting consensus per pattern
preys <- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/epi_NLR.csv", row.names = 1)
preys<- preys[preys$consensus!="",]
preys$group<- preys$consensus
chrs = unique(preys$consensus)


#blumeria
thr <- 0.003
ann = read.table("~/iowa_state/lab/RNAseq/DESeq2 analysis HvV2/bgh_ef1_ef2_CSEP.tsv", header = T)[,c(2,11)]
colnames(ann)<- c("gene","description")
#if plotting all epi effectors or all effectors
#preys <- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/epi_effectors.csv", row.names = 1)

preys<- data.frame(gene=read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/effectors.csv")[,9])
length(unique(preys$gene))
preys$group<- "effector"
chrs = unique(preys$group)

#Differential expression results
hvu_de[is.na(hvu_de)]<-1
#Put it all together
df <- subset(df, df$gene %in% grep(paste(unique(preys$gene), collapse="|"), df$gene, value = TRUE))
data = merge(ann, df, all.y = T)
data = merge(data, hvu_de, all.x = TRUE)
data[is.na(data)]<- 1

genotypes=c("wt","mla6","rar3", "bln1", "dm")

for(m in 1:length(chrs)){
  dc = subset(data, data$gene %in% preys$gene[preys$group==chrs[m]])  # select the preys that correspond to the bait m
  
  n = nrow(dc)
  print(n)
  
  filename = paste0("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/NLR_effector/mini_graphs_", paste(genotypes, collapse = "_"), "_", chrs[m], ".pdf")
  
  pdf(filename, width = 2, height = 2, fonts = "ArialMT", pointsize = 11)
  for(i in 1:n){
    
    # valerias table
    # symbs = dc[i,grep("wt", colnames(dc)[93:152], value = TRUE)]  # valerias table
    # symbs[is.na(symbs)] = 1
    # symbs[symbs > 0.001] = 1
    # symbs[symbs <= 0.001] = 8
    # 
    # #Assigning symbols to each genotype
    # symbs.bln1 = symbs[c(3,7,11,15,19,23)]
    # symbs.dm = symbs[c(4,8,12,16,20,24)]
    # symbs.mla6 = symbs[c(1,5,9,13,17,21)]
    # symbs.rar3 = symbs[c(2,6,10,14,18,2)]
    # #symbs.wt = rep(16, 6)
    # symbs.wt = rep(1, 6)
    
    #lettercodes
    symbs2 = dc[i, 93:152]
    letterCodes <- data.frame()
    for (time in paste0("t", c(0, 16, 20, 24, 32, 48))) {
      symbs2_time<- symbs2[,c(grep(time, colnames(symbs2), value = TRUE))]
      vector_multC <- as.numeric(symbs2_time[1,])
      names(vector_multC) <- apply(data.frame(colnames(symbs2_time)), 1, FUN= function(x){
        paste0(strsplit(x, "_")[[1]][2], "-", strsplit(x, "_")[[1]][3])})
      compLetters <- multcompLetters(vector_multC, threshold = thr)
      if(ncol(letterCodes)==0){
        letterCodes <- data.frame(compLetters$Letters, stringsAsFactors=FALSE)
      }else{
        letterCodes <- cbind(letterCodes, data.frame(compLetters$Letters, stringsAsFactors=FALSE))}
      colnames(letterCodes)[ncol(letterCodes)]<- time
      rm(vector_multC, compLetters)
    }
    
    
    gi = dc[i,3:92]
    
    # Calculating means for gi
    gmeans = c()
    j = 1
    for(k in 1:30){
      gmeans[k] = mean(as.numeric(gi[1,j:(j+2)]))
      j = j+3
    }
    rm(j, k, gi)
    
    # For V2 new annotation 
    means.bln1.norm = gmeans[19:24]
    means.dm.norm = gmeans[25:30]
    means.mla6.norm = gmeans[7:12]
    means.rar3.norm = gmeans[13:18]
    means.wt.norm = gmeans[1:6]
    
    # Calculating y limit
    down2=min(gmeans[1:30])
    if(is.na(down2) == TRUE){down2 = 0}
    up2=max(gmeans[1:30])
    if(is.na(up2) == TRUE){up2 = down2}
    lim2=c(down2, up2)
    
    
    # Creating a nice title
    title = dc[i,c("gene")]
    par(mgp=c(.5,.5,0), mar = c(1.5, 1.5,1.3, .1))
    # Plotting
    # plotCI(c(1,5,6,7,9,13), unlist(means.mla6.norm), lty=2, xaxt="n", xlim=c(0.5,13.5), ylim=as.numeric(lim2), gap=0, ylab="DESeq2 Normalized Read Counts by Timepoint", xlab="Hours after Inoculation", main=title, type="n", col="red")
    # lines(c(1,5,6,7,9,13), unlist(means.mla6.norm), col="red", lty=2)
    # text(c(1,5,6,7,9,13), unlist(means.mla6.norm), as.character(letterCodes["mla6",]), col="red")
    # 
    # plotCI(c(0.8,4.8,5.8,6.8,8.8,12.8), unlist(means.wt.norm), lty=1, xaxt="n", xlim=c(0.3,13.7), ylim=as.numeric(lim2), gap=0, add=TRUE, type="n", col="black")
    # lines(c(0.8,4.8,5.8,6.8,8.8,12.8), unlist(means.wt.norm), col="black", lty=1)
    # text(c(0.8,4.8,5.8,6.8,8.8,12.8), unlist(means.wt.norm), as.character(letterCodes["wt",]), col="black")
    # 
    
    
    if("wt" %in% genotypes){
      plotCI(c(0.8,4.8,5.8,6.8,8.8,12.8), unlist(means.wt.norm), lty=1, xaxt="n", xlim=c(0.3,13.7), ylim=as.numeric(lim2), gap=0, ylab="", xlab="", main=title, cex.main=0.5, type="n", col="black", lwd = 1.5)
      lines(c(0.8,4.8,5.8,6.8,8.8,12.8), unlist(means.wt.norm), col="black", lty=1)
      text(c(0.8,4.8,5.8,6.8,8.8,12.8), unlist(means.wt.norm), as.character(letterCodes["wt",]), col="black")
      
    }
    if("mla6" %in% genotypes){
      plotCI(c(1,5,6,7,9,13), unlist(means.mla6.norm), lty=2, xaxt="n", xlim=c(0.5,13.5), ylim=as.numeric(lim2), gap=0, add=TRUE, lwd = 1.5,type="n", col="red")
      lines(c(1,5,6,7,9,13), unlist(means.mla6.norm), col="red", lty=2)
      text(c(1,5,6,7,9,13), unlist(means.mla6.norm), as.character(letterCodes["mla6",]), col="red")
      
    }
    if("bln1" %in% genotypes){
      plotCI(c(0.9,4.9,5.9,6.9,8.9,12.9), unlist(means.bln1.norm), lty=1, xaxt="n", xlim=c(0.3,13.7), ylim=as.numeric(lim2), gap=0, add=TRUE, lwd = 1.5,type="n", col="brown")
      lines(c(0.9,4.9,5.9,6.9,8.9,12.9), unlist(means.bln1.norm), col="brown", lty=1)
      text(c(0.9,4.9,5.9,6.9,8.9,12.9),  unlist(means.bln1.norm), as.character(letterCodes["bln1",]), col="brown")
      
    }
    if("rar3" %in% genotypes){
      plotCI(c(1.1,5.1,6.1,7.1,9.1,13.1), unlist(means.rar3.norm), lty=2, xaxt="n", xlim=c(0.3,13.7), ylim=as.numeric(lim2), gap=0, add=TRUE, lwd = 1.5,type="n", col="darkgreen")
      lines(c(1.1,5.1,6.1,7.1,9.1,13.1), unlist(means.rar3.norm), col="darkgreen", lty=2)
      text(c(1.1,5.1,6.1,7.1,9.1,13.1), unlist(means.rar3.norm), as.character(letterCodes["rar3",]), col="darkgreen")
      
    }
    if("dm" %in% genotypes){
      plotCI(c(1.2,5.2,6.2,7.2,9.2,13.2), unlist(means.dm.norm), lty=2, xaxt="n", xlim=c(0.3,13.7), ylim=as.numeric(lim2), gap=0, add=TRUE, lwd = 1.5,type="n", col="blue")
      lines(c(1.2,5.2,6.2,7.2,9.2,13.2), unlist(means.dm.norm), col="blue", lty=2)
      text(c(1.2,5.2,6.2,7.2,9.2,13.2), unlist(means.dm.norm), as.character(letterCodes["dm",]), col="blue")
      
    }
    axis(side=1, at=c(1,5,6,7,9,13), labels=c(0,16,20,24,32,48), cex=0.5, lwd = 1.5)
    #legend("top", legend=c("CI 16151", expression(italic("bln1-m19089")), expression(italic("mla6-m18982")), expression(italic("rar3-m11526")), expression(italic("(mla6+bln1)-m19028"))),lty=c(1,1,2,2,2), col=c("black","brown","red","darkgreen","blue"), bty = "n", horiz = TRUE, text.col = c("black","brown","red","darkgreen","blue"), cex = 0.75)
    
    rm(gmeans, lim2)
  }
  
  
  dev.off()
}




#matrix of genes

library(multcompView)
hvu_de = read.csv("~/iowa_state/lab/RNAseq/DESeq2 HvR3/hv_R3_genes_deseq2_results_pairwise_genotype_padj_tax_sp.csv")
hvu_de[is.na(hvu_de)]<-1

#create a function to make the matrix with the genes and the comparisons 
make_comp_matrix<- function(de_table, gene_list, thr=0.001){
  horvu_DE_table<- de_table[grep(paste(gene_list, collapse="|"), de_table$gene), ]
  n<- nrow(horvu_DE_table)
  result<- data.frame(matrix(nrow = length(rownames(horvu_DE_table)), 
                             ncol = length(colnames(horvu_DE_table)[-1]), 
                             dimnames = list(horvu_DE_table$gene, colnames(horvu_DE_table)[-1])))
  #n=1
  for(i in 1:n){
    # I will calculate the significance groups across genotype/timepoint comparisons
    symbs2 = horvu_DE_table[i, 2:61]
    letterCodes <- c()
    for (time in paste0("t", c(0, 16, 20, 24, 32, 48))) {
      symbs2_time<- symbs2[,c(grep(time, colnames(symbs2), value = TRUE))]
      vector_multC <- as.numeric(symbs2_time[1,])
      names(vector_multC) <- apply(data.frame(colnames(symbs2_time)), 1, FUN= function(x){
        paste0(strsplit(x, "_")[[1]][2], "-", strsplit(x, "_")[[1]][3])})
      compLetters <- multcompLetters(vector_multC, threshold = thr)
      if(is.null(letterCodes)){
        letterCodes <- compLetters$Letters
      }else{letterCodes <- c(letterCodes, compLetters$Letters)}
      rm(vector_multC, compLetters)
    }
    result[i,]<-letterCodes
  }
  return(result)
}

NLR_mat<- make_comp_matrix(de_table=hvu_de, gene_list=unique(epi_NLR$gene), thr=0.001)


padj_table<- read.csv("~/iowa_state/lab/RNAseq/DESeq2 HvR3/hv_R3_genes_deseq2_results_pairwise_genotype_padj_tax_sp.csv", row.names = 1)
#padj_table<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/hv_R2_genes_deseq2_results_pairwise_genotype_padj_tax_sp.csv", row.names = 1)
padj_table_horvu<- padj_table[grep("HORVU", rownames(padj_table)), grep("wt", colnames(padj_table))]
#padj_table_horvu<- padj_table[grep("HORVU", rownames(padj_table)),]
padj_table_horvu_de<- padj_table_horvu %>% filter(if_any(everything(),`<`, 0.001))
padj_table_horvu_de<- apply(padj_table_horvu_de, 2, FUN=function(x)ifelse(x>0.001 | is.na(x), NA,1))
sum(rowSums(padj_table_horvu_de, na.rm = T)==0)
hist(rowSums(padj_table_horvu_de, na.rm = T))

padj_table_horvu_de2 <- padj_table_horvu_de[selected_genes, ]

fc_table<- read.csv("~/iowa_state/lab/RNAseq/DESeq2 HvR3/hv_R3_genes_deseq2_results_pairwise_genotype_logfc_tax_sp.csv", row.names = 1)
#fc_table<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/hv_R2_genes_deseq2_results_pairwise_genotype_logfc_tax_sp.csv", row.names = 1)
fc_table2<- fc_table[rownames(padj_table_horvu_de2), grep("wt", colnames(fc_table))]
fc_table2<- apply(fc_table2, 2, FUN=function(x)ifelse(x>0, 1,-1))
mat_horvu<-fc_table2*padj_table_horvu_de2
mat_horvu<- apply(mat_horvu, 2, FUN=function(x)ifelse(is.na(x), 0,x))
mat_horvu<- mat_horvu[order(rowSums(mat_horvu, na.rm = T), decreasing=F),]
label<- substr(colnames(padj_table_horvu_de2),1,nchar(colnames(padj_table_horvu_de2))-5)

pdf(paste0("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/transcriptome characterization/HVR3/mat_DEgenes_all2.pdf"), width = 20, height = 30, fonts = "ArialMT", pointsize = 30)
color_matrix(padj_table_horvu_de2,main="DE genes barley",show.values=F,axes=FALSE,xlab = "Timepoint/Genotype", ylab = "Gene",
             cellcolors = "red"[match(padj_table_horvu_de2, sort(unique(as.vector(padj_table_horvu_de2))))])
axis(1,at=0.5:23.5,labels=label, las=2, cex.axis = 0.5)
# color_matrix(mat_horvu,main="DE genes barley",show.values=F,axes=FALSE, border = "white",na.color = "white", xlab = "Timepoint/Genotype", ylab = "Genotype",
#              cellcolors = c("magenta", "green")[match(mat_horvu, sort(unique(as.vector(mat_horvu))))])
# 
# color_matrix(mat_horvu,main="DE genes barley",show.values=F,axes=FALSE, na.color = "black", xlab = "Timepoint/Genotype", ylab = "Genotype",
#              cellcolors = c("magenta", "green")[match(mat_horvu, sort(unique(as.vector(mat_horvu))))])
color_matrix(mat_horvu,main="DE genes barley",show.values=F,axes=FALSE, border = "black",na.color = "black", xlab = "", ylab = "Gene",
             cellcolors = c("green", "black","magenta")[match(mat_horvu, sort(unique(as.vector(mat_horvu))))])
axis(1,at=0.5:23.5,labels=label, las=2, cex.axis = 0.8)

dev.off()







