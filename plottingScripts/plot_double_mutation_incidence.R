#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)
library(binom)
library(ggrepel)


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

plot_mutation_count_histogram <- function(df, col, title, ymax){
  df <- my.rename(df, col, "xAxisVal")
  colors <- c(rep("gray",2), rep("black",6))
  plt <- ggplot(df, aes(x=xAxisVal))+
    geom_histogram(binwidth=1, fill=colors)+
    ggtitle(title)+
    ylim(0,ymax)+
    xlim(-1,6)+
    ylab('n cases')+
    xlab('N oncogenic mutations in gene')
  return(plt)
}

plot_double_mutation_comparissons <- function(df, segmentStartXParam, segmentStartYParam, segmentEndXParam, segmentEndYParam,
                                              yLabel='Fraction of cases with oncogenic multiples', xLabel='Fraction of oncogenic gene mutated cases that are oncogenic multiplets', title=''){
  df <- my.rename(df, segmentStartXParam, "segmentStartX")
  df <- my.rename(df, segmentStartYParam, "segmentStartY")
  df <- my.rename(df, segmentEndXParam, "segmentEndX")
  df <- my.rename(df, segmentEndYParam, "segmentEndY")
  
  plt <- ggplot(df, aes(x=segmentStartX, y=segmentStartY, label=textLabel))+
    geom_text(colour='red')+
    #geom_point(colour='red', size=0.75)+
    geom_text(aes(x=segmentEndX, y=segmentEndY), colour='blue')+
    #geom_point(aes(x=segmentEndX, y=segmentEndY), colour='blue', size=0.75)+
    geom_segment(aes(xend=segmentEndX, yend=segmentEndY), alpha=0.4)+
    xlim(0,1)+
    ylim(0,1)+
    xlab(xLabel)+
    ylab(yLabel)+
    ggtitle(title)
  
  return(plt)
}

plot_double_mutation_by_gene <- function(df){
  plt <- ggplot(df, aes(x=reorder(geneClassification, orderingVal), y=frac_cohort_multiplet_oncogenic_hypermutator))+
    geom_boxplot(fatten = NULL)+ #TODO FIX THIS TO PLOT THE MIDDLE LINE AS THE MEAN!!!
    stat_summary(fun.y = median, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1, linetype = "solid")+
    geom_text(aes(label=displayGeneName))+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
    ggtitle('Double mutation incidence as a function of Gene Mutation type')+
    ylab('Fraction of hypermutators of cancer type with double oncogenic mutation')+
    xlab('Gene Classification and Cancer Type')
  plt <- plt+ geom_jitter(shape=16, position=position_jitter(0.1))
  plt <- plt + labs(caption = "data:compare_double_mutation_incidence.py; plotting: plot_double_mutation_incidence.R")
  return(plt)
}

#OLD INCIDENCE METHOD

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


#################

#DIFFERENT VERSION WITH mutation graphs
df <- read.table('~/Desktop/WORK/dataForLocalPlotting/mutationNumberDistributions.tsv',sep = '\t', header=TRUE)

arrPltEndo <- plot_grid(
plot_mutation_count_histogram(df, "n_PTEN_muts", "PTEN", dim(df)[1]),
plot_mutation_count_histogram(df, "n_PIK3CA_muts", "PIK3CA", dim(df)[1]),
plot_mutation_count_histogram(df, "n_ARID1A_muts", "ARID1A", dim(df)[1]),
plot_mutation_count_histogram(df, "n_TP53_muts", "TP53", dim(df)[1]),
plot_mutation_count_histogram(df, "n_NF1_muts", "NF1", dim(df)[1]),
plot_mutation_count_histogram(df, "n_FAT1_muts", "FAT1", dim(df)[1]),
plot_mutation_count_histogram(df, "n_BRCA2_muts", "BRCA2", dim(df)[1]),
plot_mutation_count_histogram(df, "n_PIK3R1_muts", "PIK3R1", dim(df)[1]),
nrow = 2, ncol = 4
)

#arrPlt <- plot_grid(pten, pik3ca, arid1a, tp53, nf1, fat1, brca2, apc, nrow = 2, ncol = 4)

ggsave('~/Desktop/pltEnd.pdf', plot=arrPltEndo,  width = 15, height = 10, units = c("in"), limitsize = FALSE)


df <- read.table('~/Desktop/WORK/dataForLocalPlotting/mutationNumberDistributions2.tsv',sep = '\t', header=TRUE)

df$PDB

arrPltColorectal <- plot_grid(
  plot_mutation_count_histogram(df, "n_APC_muts", "APC", dim(df)[1]),
  plot_mutation_count_histogram(df, "n_PIK3CA_muts", "PIK3CA", dim(df)[1]),
  plot_mutation_count_histogram(df, "n_KRAS_muts", "KRAS", dim(df)[1]),
  plot_mutation_count_histogram(df, "n_TP53_muts", "TP53", dim(df)[1]),
  plot_mutation_count_histogram(df, "n_BRAF_muts", "BRAF", dim(df)[1]),
  plot_mutation_count_histogram(df, "n_FBXW7_muts", "FBXW7", dim(df)[1]),
  plot_mutation_count_histogram(df, "n_SMAD4_muts", "SMAD4", dim(df)[1]),
  plot_mutation_count_histogram(df, "n_PTEN_muts", "PTEN", dim(df)[1]),
  plot_mutation_count_histogram(df, "n_MSH6_muts", "MSH6", dim(df)[1]),
  plot_mutation_count_histogram(df, "n_MSH2_muts", "MSH2", dim(df)[1]),
  plot_mutation_count_histogram(df, "n_MLH1_muts", "MLH1", dim(df)[1]),
  plot_mutation_count_histogram(df, "n_NF1_muts", "NF1", dim(df)[1]),
  nrow = 3, ncol = 4
)

ggsave('~/Desktop/pltColorectal.pdf', plot=arrPltColorectal,  width = 15, height = 10, units = c("in"), limitsize = FALSE)

#GLIOMA GLIOMA GLIOMA
df <-read.table('~/Desktop/WORK/dataForLocalPlotting/mutationNumberDistributions_glioma.tsv',sep = '\t', header=TRUE)

arrPltGlioma <- plot_grid(
  plot_mutation_count_histogram(df, "n_TERT_muts", "TERT", dim(df)[1]),
  plot_mutation_count_histogram(df, "n_PTEN_muts", "PTEN", dim(df)[1]),
  plot_mutation_count_histogram(df, "n_IDH1_muts", "IDH1", dim(df)[1]),
  plot_mutation_count_histogram(df, "n_TP53_muts", "TP53", dim(df)[1]),
  plot_mutation_count_histogram(df, "n_NF1_muts", "NF1", dim(df)[1]),
  plot_mutation_count_histogram(df, "n_PIK3CA_muts", "PIK3CA", dim(df)[1]),
  plot_mutation_count_histogram(df, "n_EGFR_muts", "EGFR", dim(df)[1]),
  plot_mutation_count_histogram(df, "n_NOTCH1_muts", "NOTCH1", dim(df)[1]),
  plot_mutation_count_histogram(df, "n_MSH6_muts", "MSH6", dim(df)[1]),
  plot_mutation_count_histogram(df, "n_MSH2_muts", "MSH2", dim(df)[1]),
  plot_mutation_count_histogram(df, "n_MLH1_muts", "MLH1", dim(df)[1]),
  plot_mutation_count_histogram(df, "n_PMS2_muts", "PMS2", dim(df)[1]),
  nrow = 3, ncol = 4
)

ggsave('~/Desktop/plt2.pdf', plot=arrPltGlioma,  width = 15, height = 10, units = c("in"), limitsize = FALSE)


df <-read.table('~/Desktop/WORK/dataForLocalPlotting/mutationNumberDistributions_bladder.tsv',sep = '\t', header=TRUE)

arrPltBladder <- plot_grid(
  plot_mutation_count_histogram(df, "n_TERT_muts", "TERT", dim(df)[1]),
  plot_mutation_count_histogram(df, "n_FGFR3_muts", "FGFR3", dim(df)[1]),
  plot_mutation_count_histogram(df, "n_KDM6A_muts", "KDM6A", dim(df)[1]),
  plot_mutation_count_histogram(df, "n_TP53_muts", "TP53", dim(df)[1]),
  plot_mutation_count_histogram(df, "n_CREBBP_muts", "NF1", dim(df)[1]),
  plot_mutation_count_histogram(df, "n_PIK3CA_muts", "PIK3CA", dim(df)[1]),
  plot_mutation_count_histogram(df, "n_ARID1A_muts", "ARID1A", dim(df)[1]),
  plot_mutation_count_histogram(df, "n_FBXW7_muts", "FBXW7", dim(df)[1]),
  nrow = 2, ncol = 4
)

ggsave('~/Desktop/plt2.pdf', plot=arrPltBladder,  width = 15, height = 10, units = c("in"), limitsize = FALSE)


#############################################
######################################
############LATEST NEW VERSION


##############################
#ENDOMETRIAL
dfEndometrial <-read.table('~/Desktop/WORK/dataForLocalPlotting/endometrial_double_summary.tsv',sep = '\t', header=TRUE)
pltEndometrialNonOnc <- plot_double_mutation_comparissons(dfEndometrial, segmentStartXParam='all_multipletRatioHypermutator',
                                         segmentStartYParam='frac_cohort_multiplet_hypermutator',
                                         segmentEndXParam='all_multipletRatioNonHypermutator',
                                         segmentEndYParam='frac_cohort_multiplet_notHypermutator')
  
pltEndometrialOnc <- plot_double_mutation_comparissons(dfEndometrial, segmentStartXParam='oncogenic_multipletRatioHypermutator',
                                         segmentStartYParam='frac_cohort_multiplet_oncogenic_hypermutator',
                                         segmentEndXParam='oncogenic_multipletRatioNonHypermutator',
                                         segmentEndYParam='frac_cohort_multiplet_oncogenic_notHypermutator', title='Multiplet Oncogenic mutations in Endometrial Cancer')
  
ggsave('~/Desktop/endometrialDoubleChanges.pdf', plot=pltEndometrialOnc,  width = 8, height = 8, units = c("in"), limitsize = FALSE)

############################
#COLORECTAL
#ALERT LOADS IN THE WRONG THING FIX IT
dfColorectal <-read.table('~/Desktop/WORK/dataForLocalPlotting/geneSpecificSimulationData_Colorectal_Cancer.tsv',sep = '\t', header=TRUE)
pltColorectalOnc <- plot_double_mutation_comparissons(dfColorectal, segmentStartXParam='oncogenic_multipletRatioHypermutator',
                                                       segmentStartYParam='frac_cohort_multiplet_oncogenic_hypermutator',
                                                       segmentEndXParam='oncogenic_multipletRatioNonHypermutator',
                                                       segmentEndYParam='frac_cohort_multiplet_oncogenic_notHypermutator',
                                                      title='Multiplet mutations in Colorectal Cancer')

ggsave('~/Desktop/colorectalDoubleChanges.pdf', plot=pltColorectalOnc,  width = 11, height = 11, units = c("in"), limitsize = FALSE)

########################
#GLIOMA
dfGlioma <-read.table('~/Desktop/WORK/dataForLocalPlotting/glioma_double_summary.tsv',sep = '\t', header=TRUE)
dfGlioma$textLabel
pltGlioma <- plot_double_mutation_comparissons(dfGlioma, segmentStartXParam='all_multipletRatioHypermutator',
                                                      segmentStartYParam='frac_cohort_multiplet_hypermutator',
                                                      segmentEndXParam='all_multipletRatioNonHypermutator',
                                                      segmentEndYParam='frac_cohort_multiplet_notHypermutator',
                                                      yLabel='fraction of cases with any multiplet mutations in gene',
                                                      xLabel='fraction of gene mutated cases that are multiplets',
                                                      title='Double mutations (any type) in Glioma')

ggsave('~/Desktop/gliomaDoubleChanges.pdf', plot=pltGlioma,  width = 15, height = 15, units = c("in"), limitsize = FALSE)


##############################
df <- read.table('~/Desktop/WORK/dataForLocalPlotting/pancanDoubleSummary.tsv',sep = '\t', header=TRUE)
plt<- plot_double_mutation_by_gene(df)
ggsave('~/Desktop/plt.pdf', plot=plt,  width = 12, height = 10, units = c("in"), limitsize = FALSE)

plt <- ggplot(df, aes(x=geneClassificationPanCan, y=frac_cohort_multiplet_oncogenic_hypermutator))+
    geom_boxplot(fatten = NULL)+ #TODO FIX THIS TO PLOT THE MIDDLE LINE AS THE MEAN!!!
    stat_summary(fun.y = median, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1, linetype = "solid")+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
    geom_text(aes(label=geneAndCohortDisplay))+
    ggtitle('Double mutation incidence as a function of Gene Mutation type')+
    ylab('Fraction of hypermutators with double oncogenic mutation')+
    xlab('Gene Classification')
plt<- plt+ geom_jitter(shape=16, position=position_jitter(0.1))
plt <- plt + labs(caption = "data:compare_double_mutation_incidence.py; plotting: plot_double_mutation_incidence.R")

ggsave('~/Desktop/plt.pdf', plot=plt,  width = 6, height = 6, units = c("in"), limitsize = FALSE)


df[df$frac_cohort_multiplet_oncogenic_hypermutator > .15,]$gene






ggplot(df[df$cancerType == 'Glioma',], aes(x=geneClassificationPanCan, y=nMultiplet_Hypermutator))+
  geom_boxplot(fatten = NULL, alpha=0.5)+ #TODO FIX THIS TO PLOT THE MIDDLE LINE AS THE MEAN!!!
  stat_summary(fun.y = median, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 0.75, size = 1, linetype = "solid")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
  geom_text(aes(label=displayGlioma))+
 ggtitle('genes double mutated in TMZ glioma')+
  ylab('n cases with multiplet in gene')+
  xlab('type of gene')




#############
####PHASING

df <- read.table('~/Desktop/WORK/dataForLocalPlotting/phasingInfo.tsv',sep = '\t', header=TRUE)
plt <- ggplot(df, aes(x = cancerTypeAndGene, y= frac, fill=label))+
  geom_bar(stat='identity', position = "dodge")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))
  #ylim(0,1)

ggsave('~/Desktop/threeTypesPhasing.pdf', plot=plt,  width = 15, height = 15, units = c("in"), limitsize = FALSE)


###############DOUBLE MUTATION PERMUTATION

df <- read.table('~/Desktop/WORK/dataForLocalPlotting/mutationPermutation.tsv',sep = '\t', header=TRUE)

ggplot(df, aes(x=variable, y=value))+
  stat_summary(fun.y = mean, geom = "point", aes(ymax = ..y.., ymin = ..y..),alpha=1/2)+
  geom_point(aes(y=realIncidence), colour='red', alpha=0.9)+
  geom_segment(aes(xend=variable, y=upperBound, yend=lowerBound))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  #scale_x_discrete(name ="N mutations", labels=c(1,2,3,4,5,6,7))+
  scale_y_log10()






#########################################
######################################
##################################
###############################
##########################
#######################
##################


dfEndo <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/permutationTestDoubleMutEndometrial.tsv',sep = '\t', header=TRUE)
pltEndo <- ggplot(df, aes(x=cohort, y=p))+
  geom_point()+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  geom_hline(yintercept = 0.05, colour='blue')+
  geom_text(aes(x=6, y=-.005, label='Tumor Suppresors by percentile freq in Cancer Type'))+
  ggtitle('Endomerial Cancer')+
  theme(plot.title = element_text(hjust = 0.5))


dfColo <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/permutationTestDoubleMutColorectal.tsv',sep = '\t', header=TRUE)
pltColo <- ggplot(dfColo, aes(x=cohort, y=p))+
  geom_point()+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  geom_hline(yintercept = 0.05, colour='blue')+
  geom_text(aes(x=6, y=-.005, label='Tumor Suppresors by percentile freq in Cancer Type'))+
  ggtitle('Colorectal Cancer')+
  theme(plot.title = element_text(hjust = 0.5))

dfGlio <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/permutationTestDoubleMutGlioma.tsv',sep = '\t', header=TRUE)
pltGlio <- ggplot(dfColo, aes(x=cohort, y=p))+
  geom_point()+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  geom_hline(yintercept = 0.05, colour='blue')+
  geom_text(aes(x=7, y=-.005, label='Tumor Suppresors by percentile freq in Cancer Type'))+
  ggtitle('Gliomas')+
  theme(plot.title = element_text(hjust = 0.5))

alignedPlot <- plot_grid(ggplot()+ggtitle('Permutation Test For Double Mutation Enrichment'),
                          pltEndo, pltColo,
                         ggplot()+labs(caption='gene_multiple_hit_analysis.ipynb  plot_double_mutation_incidence.R'),
                                      nrow =4, rel_heights=c(.1,1,1,.1))

ggsave('~/Desktop/plt.pdf', plot=alignedPlot,  width = 7, height = 10, units = c("in"), limitsize = FALSE)

plot_grid(nrow = 5, ncol=1)


#
####
########
#################
#####################
########################
#####################################
#####################################
################################
########################
#################
#########
#


df <- read.table('~/Desktop/WORK/dataForLocalPlotting/pancanDoubleSummary.tsv',sep = '\t', header=TRUE)
plt<- plot_double_mutation_by_gene(df)
ggsave('~/Desktop/plt.pdf', plot=plt,  width = 12, height = 10, units = c("in"), limitsize = FALSE)

plt <- ggplot(df[df$geneClassificationPanCan == 'NonRecurrentOncogene',], aes(x=geneClassificationPanCan, y=frac_cohort_multiplet_oncogenic_hypermutator))+
  geom_boxplot(fatten = NULL)+ #TODO FIX THIS TO PLOT THE MIDDLE LINE AS THE MEAN!!!
  stat_summary(fun.y = median, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 0.75, size = 1, linetype = "solid")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
  geom_text(aes(label=geneAndCohort))+
  ggtitle('Double mutation incidence as a function of Gene Mutation type')+
  ylab('Fraction of hypermutators with double oncogenic mutation')+
  xlab('Gene Classification')
plt<- plt+ geom_jitter(shape=16, position=position_jitter(0.1))
plt <- plt + labs(caption = "data:compare_double_mutation_incidence.py; plotting: plot_double_mutation_incidence.R")

ggsave('~/Desktop/plt.pdf', plot=plt,  width = 6, height = 6, units = c("in"), limitsize = FALSE)



#
###
#########
############
##############
###################
#########################
################################
########################################
###################################
##############################
########################
#####################
#################
#############
#########
######
###
#

#this is for glioma balance stuff

emptyTheme <- theme(axis.line = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/gliomaMutBalance.tsv',sep = '\t', header=TRUE)

p <- ggplot(df)+
  stat_summary(aes(x='Unbalanced', y=unbalancedOncDoubleRate), geom='bar')+
  stat_summary(aes(x='Unbalanced', y=unbalancedOncDoubleRate), geom='errorbar', width=0, colour='black')+
  stat_summary(aes(x='Balanced', y=balancedOncDoubleRate), geom='bar')+
  stat_summary(aes(x='Balanced', y=balancedOncDoubleRate), geom='errorbar',  width=0, colour='black')+
  emptyTheme+
  ylab('N Double Oncogneic Mut/Base Pair in Region')+
  xlab('Region of genome')+
  ggtitle('Double-Onc TS Glioma')

finalP <- plot_grid(p, ggplot() + labs(caption='plot_double_mutation_incidence.R, glioma_tumor_suppressor_target'), nrow=2, rel_heights = c(1, .03))
ggsave('~/Desktop/plt.pdf', plot=finalP,  width = 2.75, height = 4, units = c("in"), limitsize = FALSE)



