#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)
library(ggrepel)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/nmutUnderSelection_endometrial.tsv',sep = '\t', header=TRUE)

plottingLevelsNeutral = c(
  "Colorectal_Normal", "Colorectal_Hypermutated",
  "Endometrial_Normal", "Endometrial_Hypermutated",
  "Glioma_Normal", "Glioma_Hypermutated", "VUS"
)
colorPalette  <- c("#3CB371", "#006400",
  "#CD5C5C", "#800000",
  "#6495ED", "#00008B", "#a9a9a9"
)

names(colorPalette)= plottingLevelsNeutral

ggplot(df[df$nOnc < 30,], aes(x = nOnc, y = nmutUnderSelection, group=displayGroup, colour=displayGroup))+
  geom_smooth(method = 'glm')+
  
  #geom_point()+
  scale_x_log10(limits= c(1,30))+
  scale_y_log10(limits= c(1,30))+
  
  ylab('Log n driver mutations under selection')+
  xlab('Log n driver mutations in case')+
  scale_colour_manual(values=colorPalette)+
  ggtitle('As the number of Oncogenic mutations increases\nthe rate of mutation selection decreases')

#method='lm'
arrPlot <- plot_grid(p1, p2, ncol=2)


df$type
