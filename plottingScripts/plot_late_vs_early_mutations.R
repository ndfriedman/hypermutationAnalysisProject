#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)
library(ggrepel)
library(ggnewscale)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

dfEndo <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/subclonalSecondHitsEndo.tsv',sep = '\t', header=TRUE)
dfColo <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/subclonalSecondHitsColo.tsv',sep = '\t', header=TRUE)

emptyTheme <- theme(axis.line = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())

pltEndo <- ggplot(dfEndo, aes(x=notSecondHitRate*1e6, y=secondHitRate*1e6))+
  geom_segment(aes(x=0, xend=20, y=0, yend=20), linetype= 2)+
  geom_point(aes(colour=Tumor_Sample_Barcode))+
  geom_text_repel(aes(label=SubclonalSecondHits, colour=Tumor_Sample_Barcode), size=2)+
  theme(legend.position="none")+
  ylab('Oncogenic Mutations/Mb (Second Hit Susceptible)')+
  xlab('Oncogenic Mutations/Mb (Second Hit Not Susceptible)')+
  emptyTheme+
  ggtitle('Endometrial Cancer')

pltColo <- ggplot(dfColo, aes(x=notSecondHitRate*1e6, y=secondHitRate*1e6))+
  geom_segment(aes(x=0, xend=15, y=0, yend=15), linetype= 2)+
  geom_point(aes(colour=Tumor_Sample_Barcode))+
  geom_text_repel(aes(label=SubclonalSecondHits, colour=Tumor_Sample_Barcode), size=2)+
  theme(legend.position="none")+
  ylab('Oncogenic Mutations/Mb (Second Hit Susceptible)')+
  xlab('Oncogenic Mutations/Mb (Second Hit Not Susceptible)')+
  emptyTheme+
  ggtitle('Colorectal Cancer')

alignedPlot <- plot_grid(pltEndo, pltColo, ncol=2)
alignedPlotWithTitle <- plot_grid(ggplot()+ggtitle('Second Hits to Tumor Suppressors in Subclones') + theme(plot.title=element_text(hjust=0.5, face='bold')),
                                  alignedPlot, ggplot()+labs(caption='late_occuring_mutations_and_double_mutation.ipynb, plot_late_vs_early_mutations.R'),
                                  rel_heights = c(.1,1,.025), nrow=3)

ggsave('~/Desktop/plot.pdf', plot=alignedPlotWithTitle,  width = 10, height = 6, units = c("in"))





####
#######
#########
################
####################
#########################
##################
############
###########
#######

#MORE complicated plot
#a bunch of columns

emptyTheme <- theme(axis.line = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())

#TODO clip the nmut column for display purposes
plot_nmut_col <- function(dfN){
  
  #NOTE EVERYTHING IN THIS PLOT IS DONE WITH COORD FLIP
  plt <- ggplot(dfN, aes(x=displayName, y=nmut))+
    
    #BACKGROUND PANEL BAR
    geom_bar(stat='identity', aes(y=max(dfN$nmut + 1), fill = Tumor_Sample_Barcode), 
             width=1.1)+ #WIDHT=1.1 allows for the gray case by case alternation
    #FILL WITH INFITIELY repeating gray and white
    scale_fill_manual(values=rep(c('gray', 'white'), length(unique(dfN$Tumor_Sample_Barcode))))+
    
    new_scale_fill()+ #RESET THE COLOR SCALE
    geom_bar(stat='identity', aes(fill=clonalStatus)) +
    scale_fill_manual(values=c('orange', 'pink'))+
    geom_text(aes(label=nmut), hjust=0, size=1)+
    coord_flip()+ 
    
    scale_x_discrete(limits = rev(levels(dfN$displayName)))+ #REVERSE ordering imposed by the coordinate flip
    ggtitle('Nmut')+
    theme(legend.position = 'none')+
    emptyTheme
  return(plt)
}

plot_onc_col <- function(dfO){
  plt <- ggplot(dfO, aes(x=displayName))+
    
    #BACKGROUND PANEL BAR
    geom_bar(stat='identity', aes(y=max(dfO$nDriverMuts + 1), fill = Tumor_Sample_Barcode), 
             width=1.1)+ #WIDHT=1.1 allows for the gray case by case alternation
    #FILL WITH INFITIELY repeating gray and white
    scale_fill_manual(values=rep(c('gray', 'white'), length(unique(dfO$Tumor_Sample_Barcode))))+
    
    
    geom_bar(stat='identity', aes(y = nDriverMuts), fill='#ADD8E6')+ #ALL ONC MUTATIONS
    geom_bar(stat='identity', aes(y = nHotspots), fill='red')+ #just the hotspots
    geom_text(aes(label=nDriverMuts, y = nDriverMuts), hjust=0, size=1)+
    coord_flip()+ 
    
    scale_x_discrete(limits = rev(levels(dfO$displayName)))+ #REVERSE ordering imposed by the coordinate flip
    ggtitle('Driver Mutations')+
    theme(legend.position = 'none')+
    theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())+ #no need for y label
    emptyTheme
  return(plt)
}

plot_hotspot_events <- function(dfJ){
  plt <- ggplot(dfJ, aes(y=0.5, x=displayName))+
    
    #BACKGROUND PANEL BAR
    geom_bar(stat='identity', aes(y=1, fill = Tumor_Sample_Barcode))+ #WIDHT=1.1 allows for the gray case by case alternation
    #FILL WITH INFITIELY repeating gray and white
    scale_fill_manual(values=rep(c('gray', 'white'), length(unique(dfJ$Tumor_Sample_Barcode))))+
    
    geom_text(aes(label=hotspotGenes))+
    theme(legend.position = 'none')+
    emptyTheme+
    theme(axis.title.y=element_blank())+
    theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())+ #no need for y label
    ggtitle('Subclonal Hotspots')+
    scale_x_discrete(limits = rev(levels(dfJ$displayName)))+ #we need to reverse to follow
    coord_flip()
  return(plt)
}

plot_denovo_double_mut_events <- function(dfD){
  plt <- ggplot(dfD, aes(y=0.5, x=displayName))+
    
    #BACKGROUND PANEL BAR
    geom_bar(stat='identity', aes(y=1, fill = Tumor_Sample_Barcode))+ #WIDHT=1.1 allows for the gray case by case alternation
    #FILL WITH INFITIELY repeating gray and white
    scale_fill_manual(values=rep(c('gray', 'white'), length(unique(dfD$Tumor_Sample_Barcode))))+
    
    geom_text(aes(label=doubleMutationNames))+
    theme(legend.position = 'none')+
    emptyTheme+
    theme(axis.title.y=element_blank())+
    theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())+ #no need for y label
    ggtitle('Genes mutated 2x Subclonally')+
    scale_x_discrete(limits = rev(levels(dfD$displayName)))+ #we need to reverse to follow
    coord_flip()
  return(plt)
}

plot_second_hit_events <- function(dfS){
  plt <- ggplot(dfS, aes(y=0.5, x=displayName))+
    
    #BACKGROUND PANEL BAR
    geom_bar(stat='identity', aes(y=1, fill = Tumor_Sample_Barcode))+ #WIDHT=1.1 allows for the gray case by case alternation
    #FILL WITH INFITIELY repeating gray and white
    scale_fill_manual(values=rep(c('gray', 'white'), length(unique(dfS$Tumor_Sample_Barcode))))+
    
    geom_text(aes(label=secondHitGenes))+
    theme(legend.position = 'none')+
    emptyTheme+
    theme(axis.title.y=element_blank())+
    theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())+ #no need for y label
    ggtitle('Tumor Suppresors With Second Subclonal Hit')+
    scale_x_discrete(limits = rev(levels(dfS$displayName)))+ #we need to reverse to follow
    coord_flip()
  return(plt)
}

plot_algined_late_event_summary <- function(dfW){
  nmutCol <- plot_nmut_col(dfW)
  oncCol <- plot_onc_col(dfW)
  dfSubclonal <- dfW[dfW$clonalStatus == 'subclonal',]
  dfSubclonal <- droplevels(dfSubclonal)
  hotspotCol <- plot_hotspot_events(dfSubclonal)
  newDoubleCol <- plot_denovo_double_mut_events(dfSubclonal)
  secondHitCol <- plot_second_hit_events(dfSubclonal)
  alignedPlot <- plot_grid(nmutCol, oncCol, hotspotCol, newDoubleCol, secondHitCol, ncol=5, rel_widths=c(1,.5, .75, .75, .75))
  return(alignedPlot)
}


dfEndo <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/subclonalSecondHitsEndo.tsv',sep = '\t', header=TRUE)

alignedPlot <- plot_algined_late_event_summary(dfEndo)
ggsave('~/Desktop/plot.pdf', alignedPlot,  width = 15, height = 15, units = c("in"))




hotspotCol <- plot_hotspot_events(dfEndo[dfEndo$clonalStatus == 'subclonal',])
