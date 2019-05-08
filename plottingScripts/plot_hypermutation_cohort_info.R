#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)
library(RColorBrewer) 
library(forcats)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

plot_cancer_type_bar <- function(df){
  pal <- colorRampPalette(brewer.pal(8,"Set1"))
  plt <- 
    df %>% 
    mutate(y = forcats::fct_reorder(cancerType, as.numeric(orderingValCType), fun = mean)) %>% 
    ggplot(aes(x=1, fill=y))+
    geom_bar(position='fill')+
    scale_fill_manual(values = sample(pal(23)))
  plt <- plt +labs(caption = "data:summarize_maf_for_hypermutation_analysis.py; plotting: plot_hypermutation_cohort_info.R")
  return(plt)
} 

plot_dominant_signature_circle <- function(df, title){
  pal <- colorRampPalette(brewer.pal(8,"Set3"))
  plt <- 
    df %>% 
    mutate(y = forcats::fct_reorder(dominantSig, as.numeric(orderingValSig), fun = mean)) %>% 
    ggplot(aes(x=1, fill=y))+
    geom_bar(position='fill', colour="transparent")+
    scale_fill_manual(values = sample(pal(10))) +
    coord_polar("y", start=0)+
    theme(axis.line = element_blank(),
          axis.ticks=element_blank(),
          axis.text=element_blank(),
          axis.title=element_blank())+
    ggtitle(title)
  return(plt)
} 

df <- read.table('~/Desktop/WORK/dataForLocalPlotting/hypermutantCohortInfo.tsv',sep = '\t', header=TRUE)

cancerTypeBar <- plot_cancer_type_bar(df)

cBar <- plot_dominant_signature_circle(df[df$cancerTypeAdj == 'Colorectal_Cancer',], title ='Colorectal Cancer')
eBar <- plot_dominant_signature_circle(df[df$cancerTypeAdj == 'Endometrial_Cancer',], title ='Endometrial Cancer')
gBar <- plot_dominant_signature_circle(df[df$cancerTypeAdj == 'Glioma',], title ='Glioma')
oBar <- plot_dominant_signature_circle(df[df$cancerTypeAdj == 'other',], title ='Other')

sigBars <- plot_grid(cBar, eBar, gBar, oBar, nrow=4)

fullPlot <- plot_grid(cancerTypeBar, sigBars, ncol=2)

ggsave('~/Desktop/plot.pdf', plot=fullPlot,  width = 10, height = 5, units = c("in"))

