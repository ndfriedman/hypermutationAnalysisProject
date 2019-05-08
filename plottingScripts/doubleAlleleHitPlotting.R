#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)
library(viridis)
library(RColorBrewer) 
library(ggrepel)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

#pal <- colorRampPalette(brewer.pal(8,"Set1"))
#scale_colour_manual(values = pal(palSize))

plot_double_mut_distribution <- function(df, title='', palSize =40){
  plt <- ggplot(df, aes(x= reorder(Gene, -orderingVal), y=AlleleCount))+
    geom_bar(aes(fill=indel), colour='white', stat='identity')+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
    #theme(legend.position="none")+
    geom_text_repel(aes(label=label),
                    nudge_x = 2.5,
                    nudge_y = 2.5,
                    segment.size  = 1,
                    segment.color = "grey50",
                    direction     = "both")+
    xlab('Gene')+
    ylab('Count')+
    ggtitle(title)
  plt <- plt +labs(caption = "data:double_mutation_vaf_analysis.py; plotting: doubleAlleleHitPlotting.R")
  return(plt)
}

endometrialDf <- read.table('~/Desktop/WORK/dataForLocalPlotting/endometrialDoubleHitPlotting.tsv',sep = '\t', header=TRUE)
endoPlot <- plot_double_mut_distribution(endometrialDf, title='Endometrial Cancer', palSize = 40)

colorectalDf <- read.table('~/Desktop/WORK/dataForLocalPlotting/colorectalDoubleHitPlotting.tsv',sep = '\t', header=TRUE)
coloPlot <- plot_double_mut_distribution(colorectalDf, title='Colorectal Cancer', palSize = 60)

gliomaDf <- read.table('~/Desktop/WORK/dataForLocalPlotting/gliomaDoubleHitPlotting.tsv',sep = '\t', header=TRUE)
gliomaPlot <- plot_double_mut_distribution(gliomaDf, title='Glioma')

ggsave('~/Desktop/plot.pdf', plot=coloPlot,  width = 6, height = 4, units = c("in"))






df <- read.table('~/Desktop/WORK/dataForLocalPlotting/endometrialDoubleMutationPlotting.tsv',sep = '\t', header=TRUE)
randoPlot <- plot_double_mut_distribution(df, title='Endometrial Cancer', palSize = 40)



