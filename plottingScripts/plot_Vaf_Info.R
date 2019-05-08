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

df <- read.table('~/Desktop/WORK/dataForLocalPlotting/colorectalVafInfo.tsv',sep = '\t', header=TRUE)
coloPlot <- plot_vaf_info_with_error_bars(df, title='Hypermutated Colorectal Cancer')
ggsave('~/Desktop/plot.pdf', plot=coloPlot,  width = 6, height = 6, units = c("in"))


df <- read.table('~/Desktop/WORK/dataForLocalPlotting/endometrialVafInfo.tsv',sep = '\t', header=TRUE)
endoPlot <- plot_vaf_info_with_error_bars(df, title='Hypermutated Endometrial Cancer')
ggsave('~/Desktop/plot.pdf', plot=endoPlot,  width = 6, height = 6, units = c("in"))

