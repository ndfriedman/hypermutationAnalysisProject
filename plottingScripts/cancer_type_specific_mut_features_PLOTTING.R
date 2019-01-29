#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(viridis)
install.packages("dplyr")

library(data.table); setDTthreads(6)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

plot_waterfall_gene_fracs <- function(df, title){
  plt <- ggplot(df, aes(x=reorder(gene, orderingVal), fill=LogGeneLength))+
    geom_bar(stat="identity", aes(y=-hypermutatorGeneFrac))+
    geom_bar(stat="identity", aes(y=nonHypermutatorGeneFrac))+
    geom_hline(aes(yintercept=0))+
    scale_fill_viridis(option="plasma")+
    coord_cartesian(ylim=c(-1, 1))+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=4))+
    ggtitle(title)+
    ylab('Bottom: hypermutator frac; Top: non-hypermutator frac')+
    xlab('Gene')+
    labs(caption = "data:cancer_type_mut_features_analysis.py; plotting: cancer_type_specific_mut_features_PLOTTING.R")
  plt <- plt + guides(fill=guide_legend(title="LOG2 IMPACT CDS Gene Length"))
  return(plt)
}


colorectalDf <- read.table('~/Desktop/WORK/dataForLocalPlotting/colorectalHypermutationComparisson.tsv',sep = '\t', header=TRUE)
pltColorectal <- plot_waterfall_gene_fracs(colorectalDf, title='Percent of Colorectal Hypermutators (n=146) vs Colorectal Non-hypermuators (n=1857) with oncogenic gene mutation')
ggsave('~/Desktop/plot.pdf', plot=plt,  width = 15, height = 10, units = c("in"))
#########################
endometrialDf <- read.table('~/Desktop/WORK/dataForLocalPlotting/endometrialHypermutationComparisson.tsv',sep = '\t', header=TRUE)
pltEndometrial <- plot_waterfall_gene_fracs(endometrialDf, title='Percent of Endometrial Hypermutators (n=41) vs Endometrial Non-hypermuators (n=706) with oncogenic gene mutation')
ggsave('~/Desktop/plot.pdf', plot=plt,  width = 15, height = 10, units = c("in"))
#########################
gliomaDf <- read.table('~/Desktop/WORK/dataForLocalPlotting/gliomaHypermutationComparisson.tsv',sep = '\t', header=TRUE)
pltGlioma <- plot_waterfall_gene_fracs(gliomaDf, title='Percent of Glioma Hypermutators (n=45) vs Glioma Non-hypermuators (n=1128) with oncogenic gene mutation')
ggsave('~/Desktop/plot.pdf', plot=plt,  width = 15, height = 10, units = c("in"))
############################
melanomaDf <- read.table('~/Desktop/WORK/dataForLocalPlotting/melanomaHypermutationComparisson.tsv',sep = '\t', header=TRUE)
pltMelanoma <- plot_waterfall_gene_fracs(melanomaDf, title='Percent of Melanoma Hypermutators (n=202) vs Melanoma Non-hypermuators (n=647) with oncogenic gene mutation')
ggsave('~/Desktop/plot.pdf', plot=plt,  width = 15, height = 10, units = c("in"))
###############################
bladderDf <- read.table('~/Desktop/WORK/dataForLocalPlotting/bladderHypermutationComparisson.tsv',sep = '\t', header=TRUE)
pltBladder <- plot_waterfall_gene_fracs(bladderDf, title='Percent of Bladder Hypermutators (n=46) vs Bladder Non-hypermuators (n=661) with oncogenic gene mutation')
ggsave('~/Desktop/plot.pdf', plot=plt,  width = 15, height = 10, units = c("in"))

plt2 <-plot_grid(
  pltGlioma, pltMelanoma, pltBladder,
  pltColorectal, pltEndometrial,
  align='hv', nrow=2, ncol=3
)

ggsave('~/Desktop/alignedPlot.pdf', plot=plt2,  width = 80, height = 40, units = c("in"), limitsize = FALSE)






