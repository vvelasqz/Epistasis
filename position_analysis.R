#position analysis

#position analysis
#folder<- "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/chr analysis/R2/"
folder<- "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/chr analysis/R3/"
#start by creating the position file from gff 
barley_gff <- read.table(file="~/iowa_state/lab/genome annotation files/tritex R3/Hv_Morex.pgsb.Jul2020.gff3",
                         sep = '\t',header = F, quote = "\"", fill=TRUE, stringsAsFactors = F, skip = 9)

# barley_gff <- read.table(file="~/iowa_state/lab/genome annotation files/tritex R2/Barley_Morex_V2_gene_annotation_PGSB.HC.gff3", 
#                       sep = '\t',header = F, quote = "\"", fill=TRUE, stringsAsFactors = F, skip = 9)

barley_gff<- barley_gff[barley_gff$V3=="gene",c(1,4,5,7,9)]
barley_gff$V9<- sapply(X = barley_gff$V9, FUN=function(x) substr(x, 4,28))
colnames(barley_gff)<- paste0("V", 1:5)

write.csv(barley_gff,paste0(folder, "Horvu_table.csv"))

# horvu_V2 <- read.table('~/iowa_state/lab/genome annotation files/mloc to horvu/Horvu_table.tsv',
#                        stringsAsFactors = F)
#plot positions
#http://bioconductor.org/packages/devel/bioc/vignettes/profileplyr/inst/doc/profileplyr.html
#https://compgenomr.github.io/book/visualizing-and-summarizing-genomic-intervals.html#visualizing-intervals-on-a-locus-of-interest

library(tidyverse)
library(GenomicRanges)
horvu <- read.csv(paste0(folder, "Horvu_table.csv"),row.names = 1, stringsAsFactors = F)
len<- horvu %>% group_by(V1) %>% summarise(lenght= max(V2, V3))
names(horvu)<- (c("chr","start","end","strand","seqname"))
Seqinfo(as.character(len$V1), seqlengths=len$lenght, isCircular=NA, genome=NA)
horvu_ranges<- makeGRangesFromDataFrame(horvu,keep.extra.columns=TRUE,ignore.strand=FALSE,
                                        seqinfo=Seqinfo(as.character(len$V1), seqlengths=len$lenght, isCircular=NA, genome=NA),
                                        seqnames.field="chr",start.field="start",end.field="end",
                                        strand.field="strand",starts.in.df.are.0based=FALSE)



bins<- data.frame(chromosome= rep(unlist(strsplit(len$V1, "r"))[-which(unlist(strsplit(len$V1, "r"))=="ch")], ceiling(len$lenght/1e6)), 
                  chr= rep(len$V1, ceiling(len$lenght/1e6)),
                  number= unlist(sapply(len$lenght, function(x) sprintf("%03d",0:floor(x/1e6)))),
                  start= unlist(sapply(len$lenght, function(x) 0:floor(x/1e6)*1e6+1)), 
                  end= unlist(sapply(len$lenght, function(x) 0:floor(x/1e6)*1e6+1e6)),
                  strand="+")

bins$number<- as.character(bins$number)
bins$bin_1Mb<- sapply(sapply(bins$number, function(x) strsplit(x, "")), function(x) paste(x, collapse = ".")) 
bins$bin_1Mb<- apply(bins, 1, function(x) paste0(x[1],".", x[7])) 
bins[sapply(1:nrow(len), function(x) sum(ceiling(len$lenght/1e6)[1:x])), "end"]<- len$lenght

bins_ranges<- makeGRangesFromDataFrame(bins,keep.extra.columns=TRUE,ignore.strand=FALSE,
                                       seqinfo=Seqinfo(as.character(len$V1), seqlengths=len$lenght, isCircular=NA, genome=NA),
                                       seqnames.field="chr",start.field="start",end.field="end",
                                       strand.field="strand",starts.in.df.are.0based=FALSE)

complete_overlap = findOverlaps(query = horvu_ranges, subject = bins_ranges, type = "any", ignore.strand=T)
complete_overlap_df = data.frame(horvu[queryHits(complete_overlap),], bins[subjectHits(complete_overlap),])
complete_overlap_df$bin_10Mb<- gsub('.{0,2}$', '', complete_overlap_df$bin_1Mb) 
complete_overlap_df$bin_100Mb<- gsub('.{0,4}$', '', complete_overlap_df$bin_1Mb) 
#colnames(complete_overlap_df)[11:13]<- c("bin_1Mb", "bin_10Mb", "bin_100Mb")
#write.csv(complete_overlap_df, paste0(folder, "Binned_genes_HR2.csv"), row.names = F)
write.csv(complete_overlap_df, paste0(folder, "Binned_genes_HR3.csv"), row.names = F)

genes_bin<-data.frame(table(complete_overlap_df$bin_1Mb))
bins<- merge(bins, genes_bin, by.x="bin_1Mb", by.y="Var1", all.x=T)
bins[is.na(bins$Freq),"Freq"]<-0
bins$bin_10Mb<- gsub('.{0,2}$', '', bins$bin_1Mb) 
bins$bin_100Mb<- gsub('.{0,4}$', '', bins$bin_1Mb) 

bins_1Mb<- bins %>% group_by(bin_1Mb) %>% summarise(chromosome=unique(chromosome), start=(min(start)), end=max(end), Freq=sum(Freq))
bins_10Mb<- bins %>% group_by(bin_10Mb) %>% summarise(chromosome=unique(chromosome), start=(min(start)), end=max(end), Freq=sum(Freq))
bins_100Mb<- bins %>% group_by(bin_100Mb) %>% summarise(chromosome=unique(chromosome), start=(min(start)), end=max(end), Freq=sum(Freq))
bins_chr<- bins %>% group_by(chromosome) %>% summarise(start=(min(start)), end=max(end), Freq=sum(Freq))

write.csv(bins_1Mb, paste0(folder, "bins_1Mb.csv"), row.names = F)
write.csv(bins_10Mb, paste0(folder, "bins_10Mb.csv"), row.names = F)
write.csv(bins_100Mb, paste0(folder, "bins_100Mb.csv"), row.names = F)
write.csv(bins_chr, paste0(folder, "bins_chrs.csv"), row.names = F)

#HG test
HG_test<- function(trts, bin_table, overlap_table){
  #table trts by bins
  hypergeom_test <- data.frame(matrix(ncol = 1+2*length(trts), nrow = nrow(bin_table))) 
  hypergeom_test[,1]<- bin_table[,1]
  colnames(hypergeom_test)<- c("binID", unlist(sapply(trts, function(x) c(paste0(x,"_pvalue"), paste0(x,"_padj")))))
  
  for(m in 1:length(trts)){
    epi_gene_list<- read.csv(paste0(folder,"Genelists for positioning/",trts[m], "_Genelist.csv"), 
                             stringsAsFactors = F, row.names = 1)[,"gene"]
    overlap_table_filtered<- overlap_table[overlap_table$seqname %in% epi_gene_list, c("seqname", colnames(bin_table)[1])]
    DE_genes_bins<-data.frame(table(overlap_table_filtered[,colnames(bin_table)[1]]))
    overlap_table_filtered<- merge(overlap_table_filtered, DE_genes_bins, by.x=colnames(bin_table)[1], by.y="Var1", all.x=T)
    DE_bins<- bin_table[unlist(bin_table[,1]) %in% DE_genes_bins$Var1,]
    
    #total DE in bin, total DE trt, sum not DE in bins listed in trt, genes in bin
    hypergeom_test[match(DE_genes_bins$Var1, unlist(hypergeom_test$binID)), 
                   paste0(trts[m],"_pvalue")]<- phyper(DE_genes_bins$Freq, length(epi_gene_list), 
                                                       sum(DE_bins$Freq - DE_genes_bins$Freq), 
                                                       DE_bins$Freq, lower.tail = FALSE)
    
    hypergeom_test[unlist(hypergeom_test$binID) %in% DE_genes_bins$Var1,paste0(trts[m],"_padj")]<- p.adjust(hypergeom_test[unlist(hypergeom_test$binID) %in% DE_genes_bins$Var1,paste0(trts[m],"_pvalue")], method = "BH")
  }
  hypergeom_test<- hypergeom_test[rowSums(is.na(hypergeom_test[,-1])) != ncol(hypergeom_test[,-1]), ]
  return(hypergeom_test)
}

#test for all bin tables
chrs = list.files(path = paste0(folder,"Genelists for positioning/"), pattern = "csv")
#trts<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/chr analysis/R2/trts_barley.csv", row.names = 1)
trts<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/chr analysis/R3/trts_barley.csv", row.names = 1)
trts<- trts$x
#trts<- gsub('.{0,13}$', '', chrs) 
#overlap_table<- read.csv(paste0(folder, "Binned_genes.csv")
#overlap_table<- read.csv(paste0(folder, "Binned_genes_HR3.csv")
complete_overlap_df<- read.csv(paste0(folder, "Binned_genes_HR3.csv"))

overlap_table<- complete_overlap_df

# bins_1Mb<- read.csv( paste0(folder, "bins_1Mb.csv"))
# bins_10Mb<-read.csv(paste0(folder, "bins_10Mb.csv"))
# bins_100Mb<-read.csv( paste0(folder, "bins_100Mb.csv"))
# bins_chr<-read.csv( paste0(folder, "bins_chrs.csv"))


HG_1Mb<- HG_test(trts, bins_1Mb, overlap_table)
write.csv(HG_1Mb, paste0(folder, "hypergeometric_test_bins_1Mb.csv"), row.names = F)

HG_10Mb<- HG_test(trts, bins_10Mb, overlap_table)
write.csv(HG_10Mb, paste0(folder, "hypergeometric_test_bins_10Mb.csv"))

HG_100Mb<- HG_test(trts, bins_100Mb, overlap_table)
write.csv(HG_100Mb, paste0(folder, "hypergeometric_test_bins_100Mb.csv"))

HG_chrs<- HG_test(trts, bins_chr, overlap_table)
write.csv(HG_chrs, paste0(folder, "hypergeometric_test_bins_chrs.csv"))


#Make the plots
#https://www.bioconductor.org/packages/release/bioc/vignettes/ggbio/inst/doc/ggbio.pdf
#BiocManager::install("ggbio")

library(ggbio)
#Now a loop
p <- autoplot(seqinfo(horvu_ranges), layout = "karyogram")
#pdf(paste0(folder, "Barley_Mla6_Bln1_chr_position_by_group.pdf"), width = 10, height = 10, fonts = "ArialMT", pointsize = 12)
pdf(paste0(folder, "Barley_chr_position_by_group.pdf"), width = 10, height = 10, fonts = "ArialMT", pointsize = 12)
for(m in 1:length(chrs)){
  name<- gsub('.{0,13}$', '', chrs[m])
  coords<- read.csv(paste0(folder,"Genelists for positioning/",chrs[m]), stringsAsFactors = F, row.names = 1)
  coords<- merge(horvu, coords, by.x="seqname",by.y="gene")
  coords_ann<-makeGRangesFromDataFrame(coords, keep.extra.columns=T, seqnames.field= "chr",start.field="start",
                                       seqinfo=Seqinfo(as.character(len$V1)), end.field="end", strand.field="strand")
  
  print(p + layout_karyogram(coords_ann, color="red") +labs(title = name))
  # if(!grepl("dmDiff", name)){
  #   print(p + layout_karyogram(coords_ann[coords_ann$Direction=="Pos",], color="red") +
  #           layout_karyogram(coords_ann[coords_ann$Direction=="Neg",], color="blue")+labs(title = name))
  # }else{
  #   print(p + layout_karyogram(coords_ann, color="purple") +labs(title = name))
  # }
}
dev.off()


#Fig 4a position
#define the significant hotspots and plot coloring by pattern
folder<- "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/chr analysis/R3/"
# HG_1Mb<- read.csv(paste0(folder, "hypergeometric_test_bins_1Mb.csv"))
# HG_10Mb<- read.csv(paste0(folder, "hypergeometric_test_bins_10Mb.csv"), row.names = 1)
HG_100Mb<- read.csv(paste0(folder, "hypergeometric_test_bins_100Mb.csv"), row.names = 1)
#HG_chrs<- read.csv(paste0(folder, "hypergeometric_test_bins_chrs.csv"), row.names = 1)

# type=c("additive_Mla6_Bln1", "symmetric", "masked", "suppression", "pseudo_masked", "positive", "negative",
#        "unique_eff_Mla6", "unique_eff_Rar3", "equal_eff_Mla6_Rar3", "predominant_Mla6", "predominant_Rar3")
type=c("additive_Mla6_Bln1", "symmetric", "masked", "suppression", "pseudo_masked", "positive", "negative")
trts<- paste0("Barley_Mla6_Bln1_",type, "_padj")

HG_100Mb_mla6_bln1<- HG_100Mb[, c(1, grep(paste(trts, collapse = "|"), colnames(HG_100Mb)))]
#thr 0.005
HG_100Mb_mla6_bln1<- HG_100Mb_mla6_bln1[rowSums(HG_100Mb_mla6_bln1[,-1] <0.005, na.rm = T)>0,
                                        c(T,colSums(HG_100Mb_mla6_bln1[,-1] <0.005, na.rm = T)>0)]
sig_trts<- sapply(colnames(HG_100Mb_mla6_bln1)[-1], FUN=function(x)strsplit(x,"_")[[1]][4])
colnames(HG_100Mb_mla6_bln1)[-1]<-sig_trts
HG_100Mb_mla6_bln1<- melt(HG_100Mb_mla6_bln1, 1, na.rm = T)
HG_100Mb_mla6_bln1<- HG_100Mb_mla6_bln1[HG_100Mb_mla6_bln1$value< 0.005,]

binned_genes<-read.csv(paste0(folder, "Binned_genes_HR3.csv"))
spread_tables_Mla6_Bln1<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/spread_tables_Mla6_Bln1.csv", row.names = 1)
spread_tables_Mla6_Bln1<- spread_tables_Mla6_Bln1[spread_tables_Mla6_Bln1$consensus %in% sig_trts, c("gene","consensus")]
spread_tables_Mla6_Bln1<- merge(spread_tables_Mla6_Bln1, binned_genes[,c(5,14)], by.x="gene", by.y="seqname", all.x=T)
sig_genes_Mla6_Bln1<- merge(spread_tables_Mla6_Bln1, HG_100Mb_mla6_bln1, by.x=c("consensus","bin_100Mb"), by.y=c("variable", "binID"))

write.csv(sig_genes_Mla6_Bln1, paste0(folder, "Barley_chr_position_significant.csv"))

library(tidyverse)
library(GenomicRanges)
horvu <- read.csv(paste0(folder, "Horvu_table.csv"),row.names = 1, stringsAsFactors = F)
len<- horvu %>% group_by(V1) %>% summarise(lenght= max(V2, V3))
names(horvu)<- (c("chr","start","end","strand","seqname"))
Seqinfo(as.character(len$V1), seqlengths=len$lenght, isCircular=NA, genome=NA)
horvu_ranges<- makeGRangesFromDataFrame(horvu,keep.extra.columns=TRUE,ignore.strand=FALSE,
                                        seqinfo=Seqinfo(as.character(len$V1), seqlengths=len$lenght, isCircular=NA, genome=NA),
                                        seqnames.field="chr",start.field="start",end.field="end",
                                        strand.field="strand",starts.in.df.are.0based=FALSE)

#Make the plots
#https://www.bioconductor.org/packages/release/bioc/vignettes/ggbio/inst/doc/ggbio.pdf
#BiocManager::install("ggbio")

library(ggbio)
#Now a loop
p <- autoplot(seqinfo(horvu_ranges), layout = "karyogram")
chrs<- sig_trts
colors  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
names(colors)<- c("additive_Mla6_Bln1", "symmetric", "masked", "suppression", "pseudo_masked", "positive", "negative")

for(m in 1:length(chrs)){
  name<- chrs[m]
  coords<- sig_genes_Mla6_Bln1[sig_genes_Mla6_Bln1$consensus==chrs[m],]
  coords<- merge(horvu, coords, by.x="seqname",by.y="gene")
  coords_ann<-makeGRangesFromDataFrame(coords, keep.extra.columns=T, seqnames.field= "chr",start.field="start",
                                       seqinfo=Seqinfo(as.character(len$V1)), end.field="end", strand.field="strand")
  p<-p + layout_karyogram(coords_ann, color=colors[chrs[m]])
}

pdf(paste0(folder, "Barley_chr_position_significant.pdf"), width = 5, height = 3, fonts = "ArialMT", pointsize = 52)
print(p)
dev.off()


#mla6 rar3
folder<- "~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/chr analysis/R3/"
HG_100Mb<- read.csv(paste0(folder, "hypergeometric_test_bins_100Mb.csv"), row.names = 1)
type=c("unique_eff_Mla6", "unique_eff_Rar3", "equal_eff_Mla6_Rar3", "predominant_Mla6", "predominant_Rar3")
trts<- paste0("Barley_Mla6_Rar3_",type, "_padj")

HG_100Mb_mla6_rar3<- HG_100Mb[, c(1, grep(paste(trts, collapse = "|"), colnames(HG_100Mb)))]
#thr 0.005
HG_100Mb_mla6_rar3<- HG_100Mb_mla6_rar3[rowSums(HG_100Mb_mla6_rar3[,-1] <0.005, na.rm = T)>0,
                                        c(T,colSums(HG_100Mb_mla6_rar3[,-1] <0.005, na.rm = T)>0)]
sig_trts<- c("unique_eff_Mla6", "equal_eff_Mla6_Rar3", "predominant_Mla6")
colnames(HG_100Mb_mla6_rar3)[-1]<-sig_trts
HG_100Mb_mla6_rar3<- melt(HG_100Mb_mla6_rar3, 1, na.rm = T)
HG_100Mb_mla6_rar3<- HG_100Mb_mla6_rar3[HG_100Mb_mla6_rar3$value< 0.005,]

binned_genes<-read.csv(paste0(folder, "Binned_genes_HR3.csv"))
spread_tables_Mla6_Rar3<- read.csv("~/iowa_state/lab/RNAseq/RNASeq_paper/vvz manuscript/analysis/epistasis_lists/barley/HVR3/spread_tables_Mla6_Rar3.csv", row.names = 1)
spread_tables_Mla6_Rar3<- spread_tables_Mla6_Rar3[spread_tables_Mla6_Rar3$consensus %in% sig_trts, c("gene","consensus")]
spread_tables_Mla6_Rar3<- merge(spread_tables_Mla6_Rar3, binned_genes[,c(5,14)], by.x="gene", by.y="seqname", all.x=T)
sig_genes_Mla6_Rar3<- merge(spread_tables_Mla6_Rar3, HG_100Mb_mla6_rar3, by.x=c("consensus","bin_100Mb"), by.y=c("variable", "binID"))

write.csv(sig_genes_Mla6_Rar3, paste0(folder, "Barley_chr_position_significant_rar3.csv"))

library(tidyverse)
library(GenomicRanges)
horvu <- read.csv(paste0(folder, "Horvu_table.csv"),row.names = 1, stringsAsFactors = F)
len<- horvu %>% group_by(V1) %>% summarise(lenght= max(V2, V3))
names(horvu)<- (c("chr","start","end","strand","seqname"))
Seqinfo(as.character(len$V1), seqlengths=len$lenght, isCircular=NA, genome=NA)
horvu_ranges<- makeGRangesFromDataFrame(horvu,keep.extra.columns=TRUE,ignore.strand=FALSE,
                                        seqinfo=Seqinfo(as.character(len$V1), seqlengths=len$lenght, isCircular=NA, genome=NA),
                                        seqnames.field="chr",start.field="start",end.field="end",
                                        strand.field="strand",starts.in.df.are.0based=FALSE)

library(ggbio)
#Now a loop
p <- autoplot(seqinfo(horvu_ranges), layout = "karyogram")
chrs<- sig_trts[-3]
colors = c("#00AFBB", "#E7B800", "#FC4E07", "#52854C", "#7570B3")
names(colors)<- c("unique_eff_Mla6", "unique_eff_Rar3", "equal_eff_Mla6_Rar3", "predominant_Mla6", "predominant_Rar3")

for(m in 1:length(chrs)){
  name<- chrs[m]
  coords<- sig_genes_Mla6_Rar3[sig_genes_Mla6_Rar3$consensus==chrs[m],]
  coords<- merge(horvu, coords, by.x="seqname",by.y="gene")
  coords_ann<-makeGRangesFromDataFrame(coords, keep.extra.columns=T, seqnames.field= "chr",start.field="start",
                                       seqinfo=Seqinfo(as.character(len$V1)), end.field="end", strand.field="strand")
  p<-p + layout_karyogram(coords_ann, color=colors[chrs[m]])
}

pdf(paste0(folder, "Barley_chr_position_significant_mla6_rar3.pdf"), width = 5, height = 3, fonts = "ArialMT", pointsize = 52)
print(p)
dev.off()


