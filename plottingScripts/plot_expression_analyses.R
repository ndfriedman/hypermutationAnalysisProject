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

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/gtexInfo.tsv',sep = '\t', header=TRUE)

p <- ggplot(df, aes(x=plottingLabel, y=rate*1e6, group=type, colour=type))+
  stat_summary(position=position_dodge(width=0.95))+
  theme(axis.text.x = element_text(angle=90))+
  ylab('Mutations/Mb')+
  xlab('Gene type & Signature Aetiology')+
  ggtitle('Effect of expression on nonsense\nmutation rate in TCGA')+
  coord_cartesian(ylim=c(0,30))+
  labs(caption = 'plot_expression_analyses.R\nanalyze_gtex_data.ipynb')

ggsave('~/Desktop/plot.pdf', plot=p,  width = 5, height = 5, units = c("in"))

#
######
############
#######
#

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/impactGenesByExpressionPole.tsv',sep = '\t', header=TRUE)
p <- ggplot(df, aes(x = plottingLabel, y=mutRate*1e6))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle=90))+
  geom_text_repel(data=df[df$mutRate*1e6 > 100,], aes(label=Gene))+
  ylab('Mutations per case per\nmegabase of gene length')+
  xlab('Gene type and expression')+
  ggtitle('Truncating mutations in POLE hypermutators\nby gene type and expression')+
  labs(caption='plot_expression_analyses.R\nanalyze_gtex_data.ipynb')

ggsave('~/Desktop/plot.pdf', plot=p,  width = 5, height = 7, units = c("in"))



