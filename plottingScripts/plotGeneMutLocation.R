#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)
require(ggridges)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

df <- read.table('~/Desktop/WORK/dataForLocalPlotting/locationDist.tsv', sep = '\t', header=TRUE)

plt <- ggplot(df, aes(y=reorder(region, -orderingVal), x=mutRate))+
  geom_density_ridges()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))

plt <- ggplot(df[df$case == 'P-0019244-T01-IM6',], aes(x=reorder(region, orderingVal), y=mutsInRegion))+
  geom_bar(stat='identity')

plt <- ggplot(df[df$case == 'P-0023123-T01-IM6',], aes(x=reorder(region, orderingVal), y=mutRate))+
  geom_bar(stat='identity')


ggsave('~/Desktop/plot.pdf', plot=plt,  width = 10, height = 10, units = c("in"))

