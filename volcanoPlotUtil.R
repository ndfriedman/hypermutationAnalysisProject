#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)
library(viridis)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

df <- read.table('~/Desktop/dataForLocalPlotting/volcanoPlotData.tsv',sep = '\t', header=TRUE)

ggplot(df, aes(x=log2oddsratio, y=log10pvalue, label=label))+
  geom_point()+
  geom_text(colour="red", size=2)+
  xlim(-2, 10)+
  ggtitle("Oncogenic/truncating gene mutations: Hypermutators vs Normal")

ggsave('~/Desktop/plot.pdf', plot=plt,  width = 1, height = 1, units = c("in"))


####TODO PUT THIS IN ITS OWN SCRIPT

barColors = c(
  "blue",
  "red",
  "green",
  "black",
  "purple",
  "orange",
  "pink",
  "yellow",
  "#2A52BE", 
  "#FF1493"
)

d2 <- read.table('~/Desktop/dataForLocalPlotting/cohortCancerTypeSummary.tsv',sep = '\t', header=TRUE)
d2$range <- factor(d2$range,levels = c("25-50mutMb", "50-75mutMb", "75-100mutMb", "100-200mutMb", "200-400mutMb", "400+mutMb"))


ggplot(d2, aes(x=range, y=value, fill=variable))+
  geom_bar(stat='identity')+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=6))+
  ylab("Number of cases")+
  scale_fill_viridis(discrete=TRUE)+
  scale_fill_manual(values=barColors)+
  ggtitle("Types of high mutation burden cancers")





