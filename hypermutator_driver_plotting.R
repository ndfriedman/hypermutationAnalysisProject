#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")


plotBoxPlot <- function(df, xValParam, yColParam, comps=NA, compPos=NA, compYPositions=NA, maxY, title, yLabelText){
  
  df <- my.rename(df, xValParam, "xVal")
  df <- my.rename(df, yColParam, "yCol")
  
  p<- ggplot(df, aes(reorder(xVal, plotOrdering), y=yCol)) +
    geom_boxplot(fatten = NULL)+ #TODO FIX THIS TO PLOT THE MIDDLE LINE AS THE MEAN!!!
    stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1, linetype = "solid")+
    coord_cartesian(ylim=c(0, maxY))+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
    ylab(yLabelText)
  
  #p<- p + stat_compare_means(comparisons = comps, label.y=compYPositions, tip.length=.01)#
  p<- p+ geom_jitter(shape=16, position=position_jitter(0.1))+
    ggtitle(title)
  
  p<- p + stat_compare_means(comparisons = 
                               comps,
                             label.y=
                               compPos, tip.length=.05)
  
  return(p)
}


df <- read.table('~/Desktop/dataForLocalPlotting/mutburdenBoxplot.tsv',sep = '\t', header=TRUE)

comparisons = 
  list( 
    c("nonApobecBladder", "APOBEC"),
    c("nonUVMelanoma", 'UV'),
    c("nonTMZGlioma", "TMZ"),
    c("ColorectalMSS", "MMR"),
    c("EndometrialMSS", "POLE"))

compPositions <- c(maxY - .05*maxY,
  maxY - .05*maxY,
  maxY - .05*maxY,
  maxY - .05*maxY,
  maxY - .05*maxY)

numberOfRows =2
numberOfColumns =3
plt<-plot_grid(
  plotBoxPlot(df, 'label', 'nHotspots', maxY=max(df$nHotspots), title='N Hotspots', yLabelText = 'n hotspot mutations'),
  plotBoxPlot(df, 'label', 'nOncogenicMutations', maxY=max(df$nOncogenicMutation) + 5, title='N Oncogenic Drivers', yLabelText = 'n oncoKb possibly oncogenic'),
  plotBoxPlot(df, 'label', 'nOncogenicOrHotspotMutations', maxY=max(df$nOncogenicOrHotspotMutations) + 5, title='N Drivers', yLabelText = 'n oncoKb possibly oncogenic or hotspot mutations'),
  plotBoxPlot(df, 'label', 'fracHotspotMutationsAtEnrichedMotif', maxY=1.05, title='Fraction hotspots at dominant sig motif', yLabelText = 'fraction hotspot mutations at Aging vs Hypermutator Sig Motif'),
  plotBoxPlot(df, 'label', 'fracOncogenicMutationsAtEnrichedMotif', maxY=1.05, title='Fraction Oncogenic mutations at dominant sig motif', yLabelText = 'fraction oncogenic mutations at Aging vs Hypermutator Sig Motif'),
  plotBoxPlot(df, 'label', 'fracDriverMutationsAtEnrichedMotif', maxY=1.05, title='Fraction Hotspot and Oncogenic mutations at dominant sig motif', yLabelText = 'fraction oncogenic+hotspot mutations at Aging vs Hypermutator Sig Motif'),
  #plotBoxPlot(df, 'label', 'Nmut', maxY=750, title='Mutation load', yLabelText = 'number of mutations in impact'),
  align='hv', nrow=numberOfRows, ncol=numberOfColumns
  )

ggsave('~/Desktop/testHypermutatorsAndDrivers.pdf', plot=plt, width = 7*numberOfColumns, height = 7*numberOfRows, units = c("in"))

df$nMutToHotspotRatio
df$Nmut

compars <- list( 
  c("TMZ", "UV"))
comparPos <- c(700)
plt2 <-plot_grid(
  plotBoxPlot(df, 'label', 'nMutToHotspotRatio', maxY=750, title='N mut to Hotspot ratio', yLabelText = 'n muts/n hotspots'),
  plotBoxPlot(df, 'label', 'nMutToOncogenicRatio', maxY=100, title='N mut to n oncogenic ratio', yLabelText = 'n muts/n oncogenic'),
  plotBoxPlot(df, 'label', 'Nmut', maxY=1000, title='number of mutations', yLabelText = 'n muts'),
  align='hv', nrow=1, ncol=3
)

ggsave('~/Desktop/testNoah.pdf', plot=plt2, width = 21, height = 7, units = c("in"))






