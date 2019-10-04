#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)
library(ggrepel)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

plot_vaf_info_with_error_bars <- function(df, title=''){
  plt <- ggplot(df, aes(x = 1 - dVafMean, y = maxVafMean))+
    geom_text_repel(aes(label=gene))+
    geom_point(aes(colour=gene))+
    geom_segment(aes(x = 1 - dVafLower, xend= 1 - dVafUpper, yend=maxVafMean, colour=gene), alpha = 0.25)+ #x error bar
    geom_segment(aes(xend = 1 - dVafMean, y=maxVafLower, yend=maxVafUpper, colour=gene), alpha = 0.25)+ #y error bar
    ylab('Ratio of Max Vaf of Gene/Case Mean Vaf')+
    xlab('1 - minimum delta VAF between two mutants')+
    ggtitle(title)+
    theme(plot.title = element_text(hjust = 0.5))
  return(plt)
}

plot_clonality_info_with_error_bars2 <- function(df, title=''){
  plt <- ggplot(df, aes(x = type, y = ccfInfo, colour=cancerType))+
    geom_segment(aes(y = lower, yend=upper, xend=type))+ #x error bar
    geom_point()+
    ylab('CCF')+
    xlab('type')+
    ggtitle(title)+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=8))+
    xlab('')
    #theme(plot.title = element_text(hjust = 0.5))
  return(plt)
}


df <- read.table('~/Desktop/WORK/dataForLocalPlotting/colorectalVafInfo.tsv',sep = '\t', header=TRUE)
coloPlot <- plot_vaf_info_with_error_bars(df, title='Hypermutated Colorectal Cancer')
ggsave('~/Desktop/plot.pdf', plot=coloPlot,  width = 6, height = 6, units = c("in"))


df <- read.table('~/Desktop/WORK/dataForLocalPlotting/endometrialVafInfo.tsv',sep = '\t', header=TRUE)
endoPlot <- plot_vaf_info_with_error_bars(df, title='Hypermutated Endometrial Cancer')
ggsave('~/Desktop/plot.pdf', plot=endoPlot,  width = 6, height = 6, units = c("in"))



################
df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/hypermutationCohortMutationClonality.tsv',sep = '\t', header=TRUE)
plt <- plot_clonality_info_with_error_bars2(df)
caption <- ggplot()+labs(caption='data=summarize_nmut_and_noncogenic_info_for_mutation_simulation.ipynb, plot=plot_Vaf_Info.R')
p2 <- plot_grid(plt, caption, nrow=2, rel_heights = c(1,.1))




df <- read.table('~/Desktop/WORK/dataForLocalPlotting/anecdoteMuts.tsv', sep = '\t', header=TRUE)

ggplot(df, aes(x=ccf, fill=label)) +
  geom_histogram(position="identity")





