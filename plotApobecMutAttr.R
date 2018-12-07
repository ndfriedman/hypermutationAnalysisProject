#written by noah friedman
#a template for R scripts for plotting

#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
library(egg)
library(dplyr)

library(cowplot)
library(data.table); setDTthreads(6)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

plot_oncogenic_muts_by_motif <- function(df){
  df[, logNmutbin := 10^(plyr::round_any(log10(Nmut + 1), 0.2)) - 1]
  df_bin <- df[, .(nOncogenicMutationsAtEnrichedMotif = mean(nOncogenicMutationsAtEnrichedMotif)), .(mutSource, logNmutbin)]
  plt<- ggplot(df_bin, aes(x=logNmutbin,y=nOncogenicMutationsAtEnrichedMotif,fill=mutSource)) +
    geom_col() +
    scale_x_log10()+
    scale_fill_manual(values=c("red", "orange", "#89CFF0", "#000080"))+
    theme(legend.text=element_text(size=3), legend.title=element_text(size=5))+
    ggtitle("Activating Mutations in Bladder Cancer")+
    xlab("N mutations (log)")+
    ylab("N oncogenic mutations")
    plt <- plt + guides(fill=guide_legend(title="Trinucleotide Signature Channel"))
  return(plt)
}

#df <- read.table('~/Desktop/dataForLocalPlotting/gliomaMutAttributionData.tsv',sep = '\t', header=TRUE)
df$nOncogenicMutationsAtEnrichedMotif

df <- fread('~/Desktop/dataForLocalPlotting/bladderMutAttributionData.tsv')
plt <- plot_oncogenic_muts_by_motif(df)

ggsave('~/Desktop/testNoahApobec.pdf',plot=plt)
