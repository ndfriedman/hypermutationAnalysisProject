#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)
library(stringr)

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/cohortDisplayFigure.tsv',sep = '\t', header=TRUE)

emptyTheme <- theme(#axis.line = element_blank(),
                    #axis.text.x = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title.x = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank()
                    )


plot_signature_layouts <- function(df, cancerType, bs, dropLegend=TRUE){
  
  #we adjust the y axis title to be aligned with dummay characters
  #maxStrLen = nchar('Cancer of Unknown Primary')
  #extraLen <- maxStrLen - nchar(cancerType)
  extraChars <- paste(replicate(45, "."), collapse = "")
  ylabel <- paste(cancerType, '\n.\n', extraChars, sep='')
  p <- ggplot(df, aes(y = reorder(Tumor_Sample_Barcode, orderingVal)))+
    geom_tile(aes(x=value, fill=Nmut_Mb))+
    theme(axis.text.y = element_blank())+
    emptyTheme+
    scale_x_discrete(breaks=bs, drop=FALSE)+
    scale_fill_viridis_c(option='magma', direction=-1, trans = "log", breaks =c(25,50,100,300))+
    ylab(ylabel)+
    theme(axis.title.y = element_text(angle = 0, hjust = 0))+
    theme(axis.text.x = element_blank())
  if(dropLegend){
    p <- p + theme(legend.position = 'none')
  }
  return(p)
}

plot_figure <- function(df){
  breakToInclude <- unique(df$value)
  plt<- plot_grid(
    plot_signature_layouts(df[df$cancerType == 'Cancer of Unknown Primary',], 'Cancer of Unknown Primary', breakToInclude),
    plot_signature_layouts(df[df$cancerType == 'Colorectal Cancer',], 'Colorectal Cancer', breakToInclude),
    plot_signature_layouts(df[df$cancerType == 'Endometrial Cancer',],'Endometrial Cancer', breakToInclude),
    plot_signature_layouts(df[df$cancerType == 'Esophagogastric Cancer',],'Esophagogastric Cancer', breakToInclude),
    plot_signature_layouts(df[df$cancerType == 'Glioma',],'Glioma', breakToInclude),
    plot_signature_layouts(df[df$cancerType == 'Head and Neck Cancer',],'Head and Neck Cancer', breakToInclude),
    plot_signature_layouts(df[df$cancerType == 'Pancreatic Cancer',],'Pancreatic Cancer', breakToInclude),
    plot_signature_layouts(df[df$cancerType == 'Prostate Cancer',],'Prostate Cancer', breakToInclude),
    plot_signature_layouts(df[df$cancerType == 'Other',],'Other', breakToInclude),
  
    plot_signature_layouts(df[df$cancerType == 'Bladder Cancer',],'Bladder Cancer', breakToInclude),
    plot_signature_layouts(df[df$cancerType == 'Lung Cancer',], 'Lung Cancer',breakToInclude),
    plot_signature_layouts(df[df$cancerType == 'Melanoma',],'Melanoma', breakToInclude),
    plot_signature_layouts(df[df$cancerType == 'Soft Tissue Sarcoma',],'Soft Tissue Sarcoma', breakToInclude),
  
                nrow=13, rel_heights = c(
                  dim(df[df$cancerType == 'Cancer of Unknown Primary',])[[1]],
                  dim(df[df$cancerType == 'Colorectal Cancer',])[[1]],
                  dim(df[df$cancerType == 'Endometrial Cancer',])[[1]],
                  dim(df[df$cancerType == 'Esophagogastric Cancer',])[[1]],
                  dim(df[df$cancerType == 'Glioma',])[[1]],
                  dim(df[df$cancerType == 'Head and Neck Cancer',])[[1]],
                  dim(df[df$cancerType == 'Pancreatic Cancer',])[[1]],
                  dim(df[df$cancerType == 'Prostate Cancer',])[[1]],
                  dim(df[df$cancerType == 'Other',])[[1]],
                  
                  dim(df[df$cancerType == 'Bladder Cancer',])[[1]],                                    
                  dim(df[df$cancerType == 'Lung Cancer',])[[1]],
                  dim(df[df$cancerType == 'Melanoma',])[[1]],
                  dim(df[df$cancerType == 'Soft Tissue Sarcoma',])[[1]]))
  return(plt)
}

plt <- plot_figure(df)
axisLabels <- ggplot(df, aes(x=value))+
  scale_x_discrete(breaks=unique(df$value), drop=FALSE)+
  emptyTheme+
  theme(axis.text.x = element_text(angle = 60, vjust = .5, size=20),
        axis.title.y = element_text(angle = 0))+
  ylab(paste(replicate(45, "."), collapse = ""))
plotWithTitleAndAxes <- plot_grid(ggplot()+ggtitle('Characteristics of Hypermutator Cohort')+theme(plot.title = element_text(hjust=.5, face='bold', size=25)),
                                  plt, axisLabels, rel_heights = c(.05,1,.1), nrow=3) 

legend <- get_legend(plot_signature_layouts(df[df$cancerType == 'Cancer of Unknown Primary',], 'Cancer of Unknown Primary', unique(df$value), dropLegend=FALSE))
plotWithLegend <- plot_grid(plotWithTitleAndAxes, legend, ncol=2, rel_widths = c(1,.25))

ggsave('~/Desktop/plot.pdf', plot=plotWithLegend,  width = 11, height = 15, units = c("in"))


"""plt1 <- ggplot(df[df$cancerType == 'Soft Tissue Sarcoma',], aes(y = reorder(Tumor_Sample_Barcode, orderingVal)))+
         geom_tile(aes(x=value, fill=Nmut_Mb))+
#facet_wrap( ~ cancerType, strip.position = "left", scales = "free_y", dir='v')+
theme(axis.text.y = element_blank())+
#theme( #FACET WRAP THEME
#   strip.background = element_blank(),
#  strip.text.x = element_blank(),)+
emptyTheme+
scale_fill_viridis_c(option='magma', direction=-1, trans = "log", breaks =c(25,50,100,300))+
ylab('')
"""
