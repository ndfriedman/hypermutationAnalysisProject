#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")


df <- read.table('~/Desktop/WORK/dataForLocalPlotting/multipleSampleClonalityPlotting.csv',sep = ',', header=TRUE)

plt <- ggplot(df, aes(x=factor(class, levels=c('shared', 'private', 'shared_oncogenic', 'private_oncogenic')), y=mean))+
  #stat_summary(fun.y = mean, geom = "point", size = 1)+
  geom_point()+
  geom_segment(aes(xend=class, y=lowerBound, yend=upperBound))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
  ggtitle('Clonality by class of mutation')+
  xlab('mutation class')+
  ylab('fraction of mutations that are clonal')+
  ylim(0,1)

plt <- plt +labs(caption = "data:hypermutationTwoSampleAnalysis.py; plotting: multipleSamplePlotting.R")

#############################
ggsave('~/Desktop/plot.pdf', plot=plt,  width = 2.5, height = 5, units = c("in"))


###PLOT the mutation types that accrue in private vs non private
df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/privateVsSharedMutType.tsv',sep = '\t', header=TRUE)
ggplot(df, aes(x= reorder(type, orderingVal), y=n))+
  stat_summary(aes(colour=timing))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
  ylab('Fraction of All Oncogenic Mutations in Class')+
  xlab('Type of oncogenic gene alteration')+
  ggtitle('Types of mutations ')




df[df$class == 'unrelated',]




