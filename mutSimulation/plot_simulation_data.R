#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

plotSimulationBarData <- function(df){
  plt <- ggplot(df, aes(x= reorder(factor(SigName), average), y=average,  fill = factor(modelName))) +
    geom_bar(stat='identity', position = "dodge2")+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=5))+
    ylab('N Oncogenic mutations per non oncogenic mutation')+
    xlab('Signature')+
    ggtitle('Simulated Acquisition of Oncogenic Mutations by Signature')
  return(plt)
}


simulationData <- read.table('~/Desktop/dataForLocalPlotting/simulationData.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util

plt<- plotSimulationBarData(simulationData)

ggsave('~/Desktop/noahDogTe.pdf', plt)