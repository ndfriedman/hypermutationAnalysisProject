#written by noah friedman
#a template for R scripts for plotting
#a script intended to compare the observed vs expected mtuation burdens of cases

library(ggplot2)
library(grid)
require(cowplot)
library(egg)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

#IDEA points, one in each color

plot_observed_and_exprected_oncogenic_mut_burden <- function(df, title){
  plt <- ggplot(df, aes(x=reorder(Tumor_Sample_Barcode, Nmut), y=value, color=variable))+
  geom_point()+
  
  geom_smooth(aes(x=reorder(Tumor_Sample_Barcode, Nmut), y=ratio, group=1), method="lm")+
  scale_y_log10()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=4))+
  theme(axis.ticks.x  = element_blank())+
  ylab("")+
  xlab("Tumor Sample Barcode ")+
  ggtitle(title)+
  scale_colour_viridis_d()
  #scale_colour_brewer(palette = "Spectral")
  
  
  return(plt)
}

df <- read.table('~/Desktop/dataForLocalPlotting/simulatedVsObservedDataTMZ.tsv',sep = '\t', header=TRUE)
p <- plot_observed_and_exprected_oncogenic_mut_burden(df, 'Number of observed oncogenic mutations vs simulated oncogenic mutations in glioma')
ggsave('~/Desktop/noahDogTe.pdf', p,  width = 20, height = 3, units = c("in"))

dfReduced <- df[df$Nmut > 30,]
p2 <- plot_observed_and_exprected_oncogenic_mut_burden(dfReduced, 'Number of observed oncogenic mutations vs simulated oncogenic mutations in high mutation burden gliomas')
ggsave('~/Desktop/noahDogTe.pdf', p2,  width = 10, height = 4, units = c("in"))



