#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
library(egg)
library(dplyr)

library(cowplot)
library(data.table); setDTthreads(6)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

plot_oncogenic_muts_by_motif <- function(df){
  df[, logNmutbin := 10^(plyr::round_any(log10(Nmut + 1), 0.2)) - 1]
  df_bin <- df[, .(nOncogenicMutationsAtEnrichedMotif = mean(nOncogenicMutationsAtEnrichedMotif)), .(mutSource, logNmutbin)]
  plt<- ggplot(df_bin, aes(x=logNmutbin,y=nOncogenicMutationsAtEnrichedMotif,fill=mutSource)) +
    geom_col() +
    scale_x_log10()+
    scale_fill_manual(values=c("red", "orange", "#89CFF0", "#000080"))+
    theme(legend.text=element_text(size=3), legend.title=element_text(size=5))+
    ggtitle("Activating Mutations in Glioma")+
    xlab("N mutations (log)")+
    ylab("N oncogenic mutations")
  plt <- plt + guides(fill=guide_legend(title="Trinucleotide Signature Channel"))
  return(plt)
}

#df <- read.table('~/Desktop/dataForLocalPlotting/gliomaMutAttributionData.tsv',sep = '\t', header=TRUE)

df <- fread('~/Desktop/dataForLocalPlotting/gliomaMutAttributionData.tsv')
plt <- plot_oncogenic_muts_by_motif(df)

ggsave('~/Desktop/testNoah.pdf',plot=plt)


####################################################################

geneMutDistDf = read.table('~/Desktop/dataForLocalPlotting/gliomaHypermutationDistribution.tsv',sep = '\t', header=TRUE)

ggplot(geneMutDistDf, aes(reorder(x=gene, -ordering),fill=mutClassification)) +
  geom_bar(stat = "count", position="stack")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=5))+
  ylab('N cases with mutation')+
  xlab('Gene (ordered by total number of mutations in gene)')+
  scale_fill_manual(values=c("#89CFF0", "orange", "#000080","red"))+
  theme(legend.text=element_text(size=3), legend.title=element_text(size=5))


########
geneMutDistDf$gene

plt <- ggplot()+
  geom_bar(data=geneMutDistDf, aes(reorder(x=gene, -ordering),fill=mutClassification), stat = "count", position="stack")+
  geom_point(data=geneMutDistDf, aes(reorder(x=gene, -ordering), y=30*ratio))+
  geom_smooth(data=geneMutDistDf, aes(reorder(x=gene, -ordering), y=30*ratio, group=1), method="lm")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=5))+
  coord_cartesian(ylim=c(0, 30))+
  ylab('N cases with mutation')+
  xlab('Gene (ordered by total number of mutations in gene)')+
  scale_fill_manual(values=c("#89CFF0", "orange", "#000080","red"))+
  theme(legend.text=element_text(size=15), legend.title=element_text(size=20))+
  theme(axis.title=element_text(size=9))+
  ggtitle("Top 50 most mutated genes in glioma hypermutators")
plt<-plt+scale_y_continuous(sec.axis = sec_axis(~ . *1/30, name = "Fraction of non-TMZ cases with gene mutated/Fraction of TMZ cases with gene mutated"))

geneMutDistDf[geneMutDistDf$ratio >1,]

p2 <- ggplot()+
  stat_summary(data=geneMutDistDf, aes(reorder(x=gene, -ordering), y=geneLength, colour="balck"), fun.y = "mean", geom = "bar")+
  theme(axis.title=element_text(size=9))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=5))+
  theme(legend.text=element_text(size=3), legend.title=element_text(size=5))
  
finalPlot <-plot_grid(plt, p2, 
                      align='hv', nrow=2,
                      rel_heights = c(1,.2)
                      )
  
ggsave('~/Desktop/testNoah1.pdf', plot=finalPlot,  width = 16, height = 12, units = c("in"))
  
  
  