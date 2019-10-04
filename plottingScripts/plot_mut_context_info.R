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



#
######
#########
#############
#################
###################
######################
###########################
#####################
#################
##############
##########
#######
#####
#

emptyTheme <- theme(axis.line = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.title.x= element_blank(), axis.ticks.x = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/mutationMotifSummary.tsv',sep = '\t', header=TRUE)

plot_cancer_type_motif_summary <- function(df, cancerType, maxY, hyperTitle){
  normalName <- paste(cancerType, '_Normal', sep="")
  hyperName <- paste(cancerType, '_Hyper', sep="")
  
  normalCohort = df[df$cohort == normalName,]
  hyperCohort = df[df$cohort == hyperName,]
  
  normalBars <- ggplot(normalCohort)+
    stat_summary(aes(x=type, y=n), geom='bar')+
    coord_cartesian(ylim=c(0, maxY))+
    theme(axis.text.x = element_text(angle=70, hjust=1))+
    emptyTheme+
    theme(legend.position = "none")+
    ggtitle('Normal')+
    ylab('N Onc Mut/Case')
  normalPanel <- make_fractional_panels(normalCohort)
  normalPlot <- plot_grid(normalBars, normalPanel, nrow=2, rel_heights=c(1,.15))
    
  hyperBars <- ggplot(hyperCohort)+
    stat_summary(aes(x=type, y=n), geom='bar')+
    coord_cartesian(ylim=c(0, maxY))+
    theme(axis.text.x = element_text(angle=70, hjust=1))+
    emptyTheme+
    ggtitle(hyperTitle)+
    theme(legend.position = "none")+
    ylab('N Onc Mut/Case')
  hyperPanel <- make_fractional_panels(hyperCohort)
  hyperPlot <- plot_grid(hyperBars, hyperPanel, nrow=2, rel_heights=c(1,.15))

  alignedPlot <- plot_grid(normalPlot, hyperPlot, legend, ncol=3,  rel_widths = c(1,1,.2))
  plotWithTitle <- plot_grid(ggplot()+ggtitle(cancerType)+theme(plot.title=element_text(hjust=.5, face='bold')),
                             alignedPlot, rel_heights = c(.05,1), nrow=2)
  return(plotWithTitle)
}

make_fractional_panels <-function(df){
  
  newDf <- data.frame(matrix(ncol = 2, nrow = 0))
  colNames <- c('type', 'fracRelated')
  names(newDf) <- colNames
  for(g in unique(df$type)){
    groupMeanRelated <- mean(df[df$type == g,]$frac, na.rm=TRUE)
    localDf <- data.frame(g, groupMeanRelated)
    names(localDf) <- colNames
    newDf <- rbind(newDf, localDf)
  }
  
  panelPlot <- ggplot(newDf)+
    geom_tile(aes(x=type, y = .5, fill=fracRelated))+
    scale_fill_gradient(low = "white", high = "black", limits=c(0,1))+
    #scale_fill_viridis_c(option='magma', direction=-1, limits = c(0,1))+
    theme(axis.line = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank(),
          legend.position = 'none'
    )+
    scale_y_continuous(breaks=c(.5), labels=c('   '))+ #dummy breaks and y axis label
    ylab('...')
  return(panelPlot)
}


endometrialPlot <- plot_cancer_type_motif_summary(df, 'Endometrial_Cancer', 25, 'POLE-hypermutated')
colorectalPlot <- plot_cancer_type_motif_summary(df, 'Colorectal_Cancer', 25, 'MSI-hypermutated')
gliomaPlot <- plot_cancer_type_motif_summary(df, 'Glioma', 25, 'TMZ-hypermutated')
alignedPlot <- plot_grid(endometrialPlot, colorectalPlot, gliomaPlot, ncol=3)
  
plotWithTitleAndCaption <- plot_grid(
  ggplot()+ggtitle('Oncogenic Mutation Characteristics in Hypermutated and Non-hypermutated Cancer')+theme(plot.title=element_text(hjust=.5, face='bold')),
          alignedPlot, ggplot()+labs(caption = 'plot_mut_context_info.R, hotspot_and_motif_analysis.ipynb'),
          rel_heights=c(.05,1,.025), nrow=3)

#set up the legend
newDf <-  data.frame(.2)
names(newDf) <- c('Frac_Cancer_Type_Related_Muts')
pForLegend <- ggplot(newDf)+
  geom_tile(aes(y=1, x=1, fill=Frac_Cancer_Type_Related_Muts))+
  scale_fill_gradient(low = "white", high = "black", limits=c(0,1))+
  guides(fill=guide_legend(title='Fraction Mutations\nIn Cancer-Type\nRelated Genes'))
  
leg <- get_legend(pForLegend)

finalPlotWithLegend <- plot_grid(plotWithTitleAndCaption, leg, ncol=2, rel_widths=c(1,.1))

ggsave('~/Desktop/plot.pdf', finalPlotWithLegend,  width = 12, height = 6, units = c("in"))



