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

#plt <- plt +labs(caption = "data:double_mutation_vaf_analysis.py; plotting: doubleAlleleHitPlotting.R")

emptyTheme <- theme(axis.line = element_blank(),
                      axis.title = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank())

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
    emptyTheme+
    xlab('Gene')+
    ylab('Count')+
    ggtitle(title)
  return(plt)
}



endometrialDf <- read.table('~/Desktop/WORK/dataForLocalPlotting/endometrialDoubleHitPlotting.tsv',sep = '\t', header=TRUE)
endoPlot <- plot_double_mut_distribution(endometrialDf, title='Endometrial Cancer', palSize = 40)

colorectalDf <- read.table('~/Desktop/WORK/dataForLocalPlotting/colorectalDoubleHitPlotting.tsv',sep = '\t', header=TRUE)
coloPlot <- plot_double_mut_distribution(colorectalDf, title='Colorectal Cancer', palSize = 60)

alignedPlot <- plot_grid(endoPlot, coloPlot, ncol=2)
alignedPlotFinal <- plot_grid(ggplot()+ggtitle('Mutations violating infinite sites hypothesis in cancer')+theme(plot.title = element_text(size = 15, face = "bold")),
                              alignedPlot,
                              ggplot()+labs(caption='mutations_violating_the_infinite_sites_model.ipynb  doubleAlleleHitPlotting.R'),
                              nrow=3, rel_heights = c(.1,1,.05))

gliomaDf <- read.table('~/Desktop/WORK/dataForLocalPlotting/gliomaDoubleHitPlotting.tsv',sep = '\t', header=TRUE)
gliomaPlot <- plot_double_mut_distribution(gliomaDf, title='Glioma')

ggsave('~/Desktop/plot.pdf', plot=alignedPlotFinal,  width = 10, height = 5, units = c("in"))

##########OTHER small plots

#TUMOR sup vs oncogene
df <- read.table('~/Desktop/WORK/dataForLocalPlotting/combinedDoubleHitPlotting.tsv',sep = '\t', header=TRUE)
p <- ggplot(df, aes(x=mutType))+
  geom_bar(stat='count')+
  emptyTheme+
  theme(plot.title = element_text(hjust = 0.5, face='bold'))+
  ggtitle('MVISH mutations')+
  ylab('N muts')+
  xlab('Gene type')+
  theme(axis.text.x = element_text(angle=70, hjust=1))

ggsave('~/Desktop/plot.pdf', plot=p,  width = 2, height = 3, units = c("in"))

##########################

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/MVISH_by_motif.tsv',sep = '\t', header=TRUE)
p <- ggplot(df, aes(x = reorder(motif, orderingVal)))+
  geom_bar(stat='count')+
  emptyTheme+
  theme(plot.title = element_text(hjust = 0.5, face='bold'))+
  ggtitle('MVISH mutations')+
  ylab('N muts')+
  xlab('Mutation Context')+
  theme(axis.text.x = element_text(angle=70, hjust=1))
ggsave('~/Desktop/plot.pdf', plot=p,  width = 2, height = 3, units = c("in"))


