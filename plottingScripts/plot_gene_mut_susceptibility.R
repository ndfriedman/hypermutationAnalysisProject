#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

layer_gene_plot <- function(df){
  plt1 <- ggplot(df, aes(x=reorder(gene, -orderingVal), y=nDistinctOncogenicSnps))+
    geom_bar(stat='identity')+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))+
    ggtitle('Number of distinct oncogenic mutations observed in genes')
  plt2 <- ggplot(df, aes(x=reorder(gene, -orderingVal), y=cdsLength))+
    geom_bar(stat='identity')+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))+
    labs(caption = "data:analyze_mutation_susceptibility.py; plotting: plot_gene_mut_susceptibility.R")
  alignedPlot <- plot_grid(
    plt1, plt2,
    align='hv', nrow=2, ncol=1
  )
  return(alignedPlot)
}

truncating_gene_plot <- function(df, title="title"){
  #plt1 <- ggplot(df, aes(x=reorder(Gene, -NPoleStopCodonPossible), y=NPoleStopCodonPossible))+
  #  geom_bar(stat='identity')+
  #  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  #  ggtitle(title)
  #plt2 <- ggplot(df, aes(x=reorder(Gene, -NPoleStopCodonPossible), y=NObsercedTruncAtMotif))+
  #  geom_bar(stat='identity')+
  #  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  #  ylab('n oncogenic mutations at POLE motif')
  #alignedPlot <- plot_grid(
  #  plt1, plt2,
  #  align='hv', nrow=2, ncol=1
  #)
  #return(alignedPlot)
  plt <- ggplot(df, aes(x=NPoleStopCodonPossible, y=NObsercedTruncAtMotif))+
    #geom_text(aes(label=Gene))+
    geom_point()+
    geom_smooth(method='lm',formula=y~x)
  return(plt)
}

#MUT susceptibility as number of distinct oncogenic mutations/size of impact region
df <- read.table('~/Desktop/WORK/dataForLocalPlotting/geneOncogenicInfo.tsv',sep = '\t', header=TRUE)
plt <- layer_gene_plot(df)

ggsave('~/Desktop/plot.pdf', plot=plt,  width = 75, height = 8, units = c("in"), limitsize = FALSE)


#Truncating mutations vs Possible stop codons
df1 <- read.table('~/Desktop/WORK/dataForLocalPlotting/truncatingInfoPole.tsv',sep = '\t', header=TRUE)

dfZ <- df1[df1$NObsercedTruncAtMotif > 0,]

p <- truncating_gene_plot(dfZ, title='Possible Stop Gain Mutations Introduced by POLE vs Observed Stop Gain Mutations Endometrial POLE (n=23)')
ggsave('~/Desktop/plot.pdf', plot=p,  width = 45, height = 8, units = c("in"), limitsize = FALSE)

ggsave('~/Desktop/plot.pdf', plot=p,  width = 45, height = 45, units = c("in"), limitsize = FALSE)



