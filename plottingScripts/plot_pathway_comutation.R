#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)
library(plyr)
library(reshape2)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

plot_cooccurence_table <- function(df, title){
  plt <- ggplot(df, aes(x=names, y=variable, fill=value)) + 
    geom_tile() +
    scale_fill_gradient(low="white") +
    geom_text( aes(label=value))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10))+
    ylab('')+
    xlab('')+
    ggtitle(title)
  return(plt)
}

df <- read.table('~/Desktop/WORK/dataForLocalPlotting/coOccurenceData.tsv',sep = '\t', header=TRUE)

plt <- plot_cooccurence_table(df, 'RTK-RAS Endometrial')
ggsave('~/Desktop/plot.pdf', plot=plt,  width = 10, height = 10, units = c("in"))

