#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

plot_volcano <- function(df){
  plt <- ggplot(df, aes(x=reorder(Gene, -Ratio), y=Ratio))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=5))+
  xlab("Gene")+
  ylab("Log2 Ratio")+
  ggtitle("Log Ratio Observed N Gene Mutations/Simulated in TMZ Gliomas")
  return(plt)
}

df <- read.table('~/Desktop/dataForLocalPlotting/tmzSimVsObservedGeneRatios.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
plt<- plot_volcano(df)
ggsave('~/Desktop/noahDogTe.pdf', plt,  width = 15, height = 5, units = c("in"))

