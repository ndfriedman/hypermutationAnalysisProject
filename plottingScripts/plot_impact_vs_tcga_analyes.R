#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

plot_impact_vs_tcga_comparison <- function(df, title){
  plt <- ggplot(df, aes(x=pVal, y=oddsRatio, colour=log(geneSize)))+
    geom_text_repel(aes(label=displayName))+
    geom_point()+
    scale_color_viridis_c(option='magma', direction=-1)+
    scale_x_log10()+
    scale_y_log10()+
    ylab('TCGA enriched        LogOddsRatio        IMPACT enriched')+
    ggtitle(title)
  return(plt)
}

dfColo <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/tcgaVsImpactGenesMutated_colon.tsv',sep = '\t', header=TRUE)
pltColo <- plot_impact_vs_tcga_comparison(dfColo, 'Colorectal Cancer')

dfEndo <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/tcgaVsImpactGenesMutated_endo.tsv',sep = '\t', header=TRUE)
pltEndo <- plot_impact_vs_tcga_comparison(dfEndo, 'Endometrial Cancer')

alignedPlot <- plot_grid(pltEndo, pltColo, ncol=2)
fullPlot <- plot_grid(ggplot()+ggtitle('TCGA vs IMPACT gene mutation comparissons') + theme(plot.title = element_text(hjust=0.5, face='bold')), 
                      alignedPlot,
                      ggplot() + labs(caption='plot_impact_vs_tcga_analyses.R, compare_impact_to_tcga_hypermutants.ipynb'),
                      nrow=3, rel_heights = c(.1,1,.05))

ggsave('~/Desktop/plot.pdf', plot=fullPlot,  width = 9, height = 5, units = c("in"))

