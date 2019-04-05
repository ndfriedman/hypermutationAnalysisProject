#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

plotSimulationBarData <- function(df){
  
  return(plt)
}


simulationData <- read.table('~/Desktop/WORK/dataForLocalPlotting/simulationData.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util

plt<- plotSimulationBarData(simulationData)

ggsave('~/Desktop/noahDogTe.pdf', plt)