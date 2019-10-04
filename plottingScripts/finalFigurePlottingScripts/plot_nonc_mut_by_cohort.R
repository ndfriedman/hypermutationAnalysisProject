#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

emptyTheme <- theme(axis.line = element_blank(),
                    #axis.text.x = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title.x = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())

plot_data <- function(df){
  plt <- ggplot(df, aes(x=reorder(cohort, orderingVal), y=nOncMuts))+
    geom_boxplot(fatten = NULL, outlier.shape=NA)+
    stat_summary(fun.y = median, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1, linetype = "solid")+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
    ggtitle('')+
    ylab('')+
    xlab('')+
    scale_colour_manual(values =  c('black', "#267574", 'gray', '#ADFF2F', "#9acd32", '#2A52BE'), name="Dominant\nSignature")+
    ylab('N Driver Mutations')+
    xlab('Cancer Type')+
    ggtitle('Driver Mutations')+
    emptyTheme
  plt <- plt+ geom_jitter(shape=16, position=position_jitter(0.1), aes(colour=dominantSignatureAdjusted), alpha=0.75)
  return(plt)
}

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/nOncMutByCohort.tsv',sep = '\t', header=TRUE)
plt <- plot_data(df)

ggsave('~/Desktop/plot.pdf', plot=plt,  width = 5, height = 5, units = c("in"))

###########

















##################


