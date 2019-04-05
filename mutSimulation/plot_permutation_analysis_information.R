#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

df <- read.table('~/Desktop/dataForLocalPlotting/',sep = '\t', header=TRUE)

ggsave('~/Desktop/plot.pdf', plot=plt,  width = 1, height = 1, units = c("in"))

