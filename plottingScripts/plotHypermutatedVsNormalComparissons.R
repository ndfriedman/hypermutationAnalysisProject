#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

plot_box_plt_with_panels <- function(df, title='', yLabel = TRUE){
  plt <- ggplot(df, aes(y=Nmut_case, x=reorder(class, -orderingVal)))+
    geom_bar(stat='identity')+
    geom_bar(aes(y=NUniqMut_case), stat='identity', fill='#D3D3D3')+
    coord_flip()+
    xlab('')+ #has to be backwards
    ggtitle(title)
  if(yLabel == FALSE){
    plt <- plt +
      theme(axis.ticks.y  = element_blank())+
      theme(axis.line.y  = element_blank())+
      theme(axis.title.y  = element_blank())+
      theme(axis.text.y  = element_blank())+
      scale_y_continuous(breaks = (c(0, 1,2,3,4,5,10,20)))+
      labs(caption = "data:compare_hypermutator_vs_normal_gene_mut_characteristics.py; plotting: plotHypermutatedVsNormalComparisson.R")
    
  }
  else{
    plt <- plt + theme(axis.title.x  = element_blank()) #also have no x label cause its redundant
  }
  return(plt)
}


df <- read.table('~/Desktop/WORK/dataForLocalPlotting/hyperVsNormalMutCharComp.tsv',sep = '\t', header=TRUE)

pltLeft <- plot_box_plt_with_panels(df[df$Hypermutated =='False',], title='Normal')
pltRight <- plot_box_plt_with_panels(df[df$Hypermutated =='True',], yLabel = FALSE, title='Hypermutated')

alignedPlot <- plot_grid(pltLeft, pltRight, ncol=2, align= 'h')
title <- ggdraw() + draw_label("Oncogenic Mutations Across Hypermutated and Non-Hypermutated Cancer", fontface='bold')

alignedPlotWithTitles <- plot_grid(title, alignedPlot, nrow=2, rel_heights = c(.1, 1))
ggsave('~/Desktop/plot.pdf', plot=alignedPlotWithTitles,  width = 8, height = 5, units = c("in"))

##########################
ggsave('~/Desktop/plot.pdf', plot=alignedPlotWithTitles,  width = 10, height = 5, units = c("in"))


