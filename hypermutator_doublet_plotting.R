#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)
library(ggpubr)


if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

plotBoxPlot <- function(df, xValParam, yColParam, comps=NA, compPos=NA, compYPositions=NA, maxY, title, yLabelText){
  
  df <- my.rename(df, xValParam, "xVal")
  df <- my.rename(df, yColParam, "yCol")
  
  p<- ggplot(df, aes(xVal, y=yCol)) +
    geom_boxplot(fatten = NULL)+ #TODO FIX THIS TO PLOT THE MIDDLE LINE AS THE MEAN!!!
    stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1, linetype = "solid")+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
    ylab(yLabelText)
  
  p<- p + stat_compare_means(comparisons = comps, label.y=compYPositions, tip.length=.01)#
  p<- p+ geom_jitter(shape=16, position=position_jitter(0.1))+
    ggtitle(title)
  
  p<- p + stat_compare_means(comparisons = 
                               comps,
                             label.y=
                               compPos, tip.length=.05)
  
  return(p)
}

#a little function that allow us to label bars with n occurences
n_fun <- function(x){
  adjust <- 0.05
  return(data.frame(y = mean(x) + adjust, label = paste0("n = ",length(x))))
}


compareCaseFractions <- function(df, xValParam, yColParam, xlabel='', ylabel='', title=''){
  df <- my.rename(df, xValParam, "xVal")
  df <- my.rename(df, yColParam, "yCol")
  plt<-ggplot(df, aes(x=xVal, y=yCol))+
    stat_summary(fun.y = mean, geom = "bar")+
    stat_summary(fun.data = n_fun, geom = "text")+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
    ggtitle(title)+
    xlab(xlabel)+
    ylab(ylabel)
  return(plt)
}

barPlot <- function(df, xValParam, colors, xlabel ='', ylabel = '', title=''){
  df <- my.rename(df, xValParam, "xVal")
  plt <- ggplot(df, aes(x=xVal, fill=factor(hotspotMutated, levels = 
                                              c('no doublets', 'other', 'multiple', 'ERBB2', 'PTEN', 'TP53', 'PIK3CA', 'multiple(inc:PIK3CA)'))))+
    geom_bar()+
    scale_fill_manual(values=colors)+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
    theme(legend.title=element_text(size=10))+
    ylab(ylabel)+
    xlab(xlabel)+
    ggtitle(title)
    #coord_cartesian(ylim=c(0, 200))
  plt <- plt + guides(fill=guide_legend(title="Type of doublet mutation"))
  return(plt)
}

#TODO PLOT BASED ON SIGNATURE MOTIF

dfDoublet <- read.table('~/Desktop/dataForLocalPlotting/hypermutatorDoubletAnalysis.tsv', sep = '\t', header=TRUE)

fractionsByCancerType <- compareCaseFractions(dfDoublet, xValParam='cancerType', yColParam='doubleHotspotPresent', ylabel='fraction of cases', title='Fraction of cases with double hotspot muts: HYPERMUTATORS (Nmut IMPACT > 75)')
mutBurdenCancerType <- compareCaseFractions(dfDoublet, xValParam='cancerType', yColParam='Nmut_Mb', xlabel='Cancer type', ylabel='Nmut_Mb')
pltCancerType <- plot_grid(fractionsByCancerType, mutBurdenCancerType, nrow=2, ncol=1, align='hv', rel_heights = c(1,.4))
ggsave('~/Desktop/plotCancerType.pdf', plot=pltCancerType,  width = 15, height = 12, units = c("in"))

fractionsBySignature <- compareCaseFractions(dfDoublet, xValParam='signatureType', yColParam='doubleHotspotPresent', ylabel='fraction of cases', title='Fraction of cases with double hotspot muts: HYPERMUTATORS (Nmut IMPACT > 75)')
mutBurdenBySignature <- compareCaseFractions(dfDoublet, xValParam='signatureType', yColParam='Nmut_Mb', xlabel='Cancer type', ylabel='Nmut_Mb')
pltSignatureType <- plot_grid(fractionsBySignature, mutBurdenBySignature, nrow=2, ncol=1, align='hv', rel_heights = c(1,.4))
ggsave('~/Desktop/plotSignatureType.pdf', plot=pltSignatureType,  width = 15, height = 12, units = c("in"))

barplotColors = c("#F5F5F5", "#FFFFCC", "#339933", "#3300FF", "#6600CC", "#FF3333", '#800000')
pltDoubletGenesCancerType  <- barPlot(dfDoublet, xValParam='cancerType', colors = barplotColors, xlabel ='Cancer Type',  ylabel = 'N cases', title='Hotspot doublets across hypermutated (Nmut/mb > 40) tumors')
pltDoubletGenesSignatureType  <- barPlot(dfDoublet, xValParam='signatureType', colors = barplotColors, xlabel ='signatureType',  ylabel = 'N cases', title='Hotspot doublets across hypermutated (Nmut/mb > 40) tumors')

#######################
pltDoubletGenesCancerType  <- barPlot(dfDoublet, xValParam='cancerSubtype', colors = barplotColors, xlabel ='Cancer Type',  ylabel = 'N cases', title='Hotspot doublets across hypermutated (Nmut/mb > 40) tumors')

####################
#######################
pltDoubletGenesCancerTypeWithNonHmutators  <- barPlot(dfDoublet, xValParam='cancerAndMutationSubtype', colors = barplotColors, xlabel ='Cancer Type',  ylabel = 'N cases', title='Hotspot doublets across hypermutated vs normal tumors.Clipped')
fractionsByCancerTypeAcrossAll <- compareCaseFractions(dfDoublet, xValParam='cancerAndMutationSubtype', yColParam='doubleHotspotPresent', ylabel='fraction of cases', title='Fraction of cases with double hotspot muts: HYPERMUTATORS vs Normal')


ggsave('~/Desktop/plotTest.pdf', plot=pltDoubletGenesCancerType,  width = 15, height = 8, units = c("in"))


