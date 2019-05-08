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

plot_trinuc_hotspot_dist <- function(df, title){
  plt <- ggplot(df, aes(x = reorder(quadNuc, quadNucCohortFrac), fill=geneType))+
    geom_bar(stat='count')+
    coord_flip()+
    geom_line(aes(y=quadNucCohortFrac*dim(df)[1], group=1), colour='black', size=2, alpha=0.25)+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=8),
          axis.text.y = element_text(size=5))+
    xlab('Trinucleotide Motif')+
    theme(legend.position = 'none')+
    ggtitle(title)
  return(plt)
}
#todo fix normalization 
dfEndometrial <- read.table('~/Desktop/WORK/dataForLocalPlotting/hotspotMutsByTrinuc_endometrial.tsv',sep = '\t', header=TRUE)
pltEndo <- plot_trinuc_hotspot_dist(df, title = 'Hypermutated Endometrial')

dfColorectal <- read.table('~/Desktop/WORK/dataForLocalPlotting/hotspotMutsByTrinuc_colorectal.tsv',sep = '\t', header=TRUE)
pltColorectal <- plot_trinuc_hotspot_dist(dfColorectal, title = 'Hypermutated Colorectal')

dfGlioma <- read.table('~/Desktop/WORK/dataForLocalPlotting/hotspotMutsByTrinuc_glioma.tsv',sep = '\t', header=TRUE)
pltGlioma <- plot_trinuc_hotspot_dist(dfGlioma, title = 'Hypermutated Glioma')

#now align everything and ass a title
title <- ggplot()+
  ggtitle('Hotspot mutations in hypermutated cancers')

caption <- ggplot()+
  labs(caption='data: compare_mutations_by_hypermutator_motif.py plot:plot_mut_context_info.R')

gridPlots <- plot_grid(pltEndo, pltColorectal, pltGlioma, ncol=3)
finalPlot <- plot_grid(title, gridPlots, caption, nrow=3, rel_heights = c(.1, 1, .05))

ggsave('~/Desktop/plot.pdf', plot=finalPlot,  width = 10, height = 5, units = c("in"))

