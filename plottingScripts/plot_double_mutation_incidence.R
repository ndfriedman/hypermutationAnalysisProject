#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)
library(binom)


if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

plot_averages_with_bar <- function(df, title){
  plt <- ggplot(df, aes(x=reorder(label, orderingVal), colour=isHypermutation))+
    geom_point(aes(y=value))+
    geom_segment(aes(xend=reorder(label, orderingVal), y=upperConf, yend=lowerConf))+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))+
    ylim(0,1.1)+
    xlab('Gene & Population')+
    ylab('Fraction of cohort')+
    ggtitle(title)+
    theme(legend.position="none")+
    labs(caption = "data:compare_double_mutation_incidence.py; plotting: plot_double_mutation_incidence.R")
  return(plt)
}


##########ENDOMETRIAL CANCER###############
dfEndometrialDouble <- read.table('~/Desktop/WORK/dataForLocalPlotting/endometrial_double_mutation_incidence.tsv',sep = '\t', header=TRUE)
dfEndometrialOneOncogenicDouble <- read.table('~/Desktop/WORK/dataForLocalPlotting/endometrial_double_with_oncogenic_mutation_incidence.tsv',sep = '\t', header=TRUE)
dfEndometrialTwoOncogenicDouble <- read.table('~/Desktop/WORK/dataForLocalPlotting/endometrial_double_with_both_oncogenic_mutation_incidence.tsv',sep = '\t', header=TRUE)

pltEndometrialWithDouble <- plot_averages_with_bar(dfEndometrialDouble, title='Double Mutation Prevalence in Endometrial Cancer')
pltEndometrialOneOncogenicDouble <- plot_averages_with_bar(dfEndometrialOneOncogenicDouble, title='Oncogenic + Other Mutation Prevalence in Endometrial Cancer')
pltEndometrialTwoOncogenicDouble <- plot_averages_with_bar(dfEndometrialTwoOncogenicDouble, title='Double Oncogenic Mutation Prevalence in Endometrial Cancer')

arrangedPlotEndometrial <- plot_grid(
  pltEndometrialWithDouble, pltEndometrialOneOncogenicDouble, pltEndometrialTwoOncogenicDouble,
  align='hv', nrow=1, ncol=3
)
ggsave('~/Desktop/plotEndometrial.pdf', plot=arrangedPlotEndometrial,  width = 20, height = 10, units = c("in"))

#################COLORECTAL CANCER########################

dfColorectalDouble <- read.table('~/Desktop/WORK/dataForLocalPlotting/colorectal_double_mutation_incidence.tsv',sep = '\t', header=TRUE)
dfColorectalOneOncogenicDouble <- read.table('~/Desktop/WORK/dataForLocalPlotting/colorectal_double_with_oncogenic_mutation_incidence.tsv',sep = '\t', header=TRUE)
dfColorectalTwoOncogenicDouble <- read.table('~/Desktop/WORK/dataForLocalPlotting/colorectal_double_with_both_oncogenic_mutation_incidence.tsv',sep = '\t', header=TRUE)

pltColorectalWithDouble <- plot_averages_with_bar(dfColorectalDouble, title='Double Mutation Prevalence in Colorectal Cancer')
pltColorectalOneOncogenicDouble <- plot_averages_with_bar(dfColorectalOneOncogenicDouble, title='Oncogenic + Other Mutation Prevalence in Colorectal Cancer')
pltColorectalTwoOncogenicDouble <- plot_averages_with_bar(dfColorectalTwoOncogenicDouble, title='Double Oncogenic Mutation Prevalence in Colorectal Cancer')

arrangedPlotColorectal <- plot_grid(
  pltColorectalWithDouble, pltColorectalOneOncogenicDouble, pltColorectalTwoOncogenicDouble,
  align='hv', nrow=1, ncol=3
)
ggsave('~/Desktop/plotColorectal.pdf', plot=arrangedPlotColorectal,  width = 20, height = 10, units = c("in"))

###############GLIOMA##############


dfGliomaDouble <- read.table('~/Desktop/WORK/dataForLocalPlotting/glioma_double_mutation_incidence.tsv',sep = '\t', header=TRUE)
dfGliomaOneOncogenicDouble <- read.table('~/Desktop/WORK/dataForLocalPlotting/glioma_double_with_oncogenic_mutation_incidence.tsv',sep = '\t', header=TRUE)
dfGliomaTwoOncogenicDouble <- read.table('~/Desktop/WORK/dataForLocalPlotting/glioma_double_with_both_oncogenic_mutation_incidence.tsv',sep = '\t', header=TRUE)

pltGliomaWithDouble <- plot_averages_with_bar(dfGliomaDouble, title='Double Mutation Prevalence in Glioma')
pltGliomaOneOncogenicDouble <- plot_averages_with_bar(dfGliomaOneOncogenicDouble, title='Oncogenic + Other Mutation Prevalence in Glioma')
pltGliomaTwoOncogenicDouble <- plot_averages_with_bar(dfGliomaTwoOncogenicDouble, title='Double Oncogenic Mutation Prevalence in Glioma')

arrangedPlotGlioma <- plot_grid(
  pltGliomaWithDouble, pltGliomaOneOncogenicDouble, pltGliomaTwoOncogenicDouble,
  align='hv', nrow=1, ncol=3
)
ggsave('~/Desktop/plotGlioma.pdf', plot=arrangedPlotGlioma,  width = 20, height = 10, units = c("in"))


##############BLADDER##############

dfBladderDouble <- read.table('~/Desktop/WORK/dataForLocalPlotting/bladder_double_mutation_incidence.tsv',sep = '\t', header=TRUE)
dfBladderOneOncogenicDouble <- read.table('~/Desktop/WORK/dataForLocalPlotting/bladder_double_with_oncogenic_mutation_incidence.tsv',sep = '\t', header=TRUE)
dfBladderTwoOncogenicDouble <- read.table('~/Desktop/WORK/dataForLocalPlotting/bladder_double_with_both_oncogenic_mutation_incidence.tsv',sep = '\t', header=TRUE)

pltBladderWithDouble <- plot_averages_with_bar(dfBladderDouble, title='Double Mutation Prevalence in Bladder Cancer')
pltBladderOneOncogenicDouble <- plot_averages_with_bar(dfBladderOneOncogenicDouble, title='Oncogenic + Other Mutation Prevalence in Bladder Cancer')
pltBladderTwoOncogenicDouble <- plot_averages_with_bar(dfBladderTwoOncogenicDouble, title='Double Oncogenic Mutation Prevalence in Bladder Cancer')

arrangedPlotBladder <- plot_grid(
  pltBladderWithDouble, pltBladderOneOncogenicDouble, pltBladderTwoOncogenicDouble,
  align='hv', nrow=1, ncol=3
)
ggsave('~/Desktop/plotBladder.pdf', plot=arrangedPlotBladder,  width = 20, height = 10, units = c("in"))


#############MELANOMA

dfMelanomaDouble <- read.table('~/Desktop/WORK/dataForLocalPlotting/melanoma_double_mutation_incidence.tsv',sep = '\t', header=TRUE)
dfMelanomaOneOncogenicDouble <- read.table('~/Desktop/WORK/dataForLocalPlotting/melanoma_double_with_oncogenic_mutation_incidence.tsv',sep = '\t', header=TRUE)
dfMelanomaTwoOncogenicDouble <- read.table('~/Desktop/WORK/dataForLocalPlotting/melanoma_double_with_both_oncogenic_mutation_incidence.tsv',sep = '\t', header=TRUE)

pltMelanomaWithDouble <- plot_averages_with_bar(dfMelanomaDouble, title='Double Mutation Prevalence in Melanoma')
pltMelanomaOneOncogenicDouble <- plot_averages_with_bar(dfMelanomaOneOncogenicDouble, title='Oncogenic + Other Mutation Prevalence in Melanoma')
pltMelanomaTwoOncogenicDouble <- plot_averages_with_bar(dfMelanomaTwoOncogenicDouble, title='Double Oncogenic Mutation Prevalence in Melanoma')

arrangedPlotMelanoma <- plot_grid(
  pltMelanomaWithDouble, pltMelanomaOneOncogenicDouble, pltMelanomaTwoOncogenicDouble,
  align='hv', nrow=1, ncol=3
)

ggsave('~/Desktop/plotMelanoma.pdf', plot=arrangedPlotMelanoma,  width = 20, height = 10, units = c("in"))


massiveAlignedPlot <- plot_grid(
  pltEndometrialWithDouble, pltEndometrialOneOncogenicDouble, pltEndometrialTwoOncogenicDouble,
  pltColorectalWithDouble, pltColorectalOneOncogenicDouble, pltColorectalTwoOncogenicDouble,
  pltGliomaWithDouble, pltGliomaOneOncogenicDouble, pltGliomaTwoOncogenicDouble,
  pltBladderWithDouble, pltBladderOneOncogenicDouble, pltBladderTwoOncogenicDouble,
  pltMelanomaWithDouble, pltMelanomaOneOncogenicDouble, pltMelanomaTwoOncogenicDouble,
  align='hv', nrow=5, ncol=3
)

ggsave('~/Desktop/massiveAligned.pdf', plot=massiveAlignedPlot,  width = 20, height = 50, units = c("in"), limitsize = FALSE)

