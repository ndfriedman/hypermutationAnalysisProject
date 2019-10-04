#written by noah friedman
#a template for R scripts for plotting
#a script intended to compare the observed vs expected mtuation burdens of cases

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(viridis)
library(ggrepel)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

#IDEA points, one in each color

plot_simple_model_observed_vs_expected <- function(df, title){
  p <- ggplot(df)+
  geom_smooth(aes(x = Nmut_x, y=Observed.Oncogenic), method=lm, formula = y ~ splines::bs(x, 3), colour="blue")+
  geom_point(aes(x = Nmut_x, y=Observed.Oncogenic), alpha = 0.1, colour="blue")+
  geom_segment(aes(x=0, y=0, xend=max(df$Nmut_x), yend=.26*max(df$Nmut_x)), colour='red')+
  #scale_x_log10()+
  #scale_y_log10()+
  ggtitle('Observed vs Expected (simple expected model)')
  return(p)
}

plot_observed_expected_curves <- function(df, title){
  plt <- ggplot(df)+
    geom_smooth(aes(x = Nmut, y=Expected.Oncogenic), method='lm', formula = y ~ x, colour="red")+
    geom_smooth(aes(x = Nmut, y=Observed.Oncogenic), method='lm', formula = y ~ x, colour="blue")+
    geom_point(aes(x = Nmut, y=Expected.Oncogenic), alpha = 0.1, colour="red")+
    geom_point(aes(x = Nmut, y=Observed.Oncogenic), alpha = 0.1, colour="blue")+
    scale_x_log10()+
    scale_y_log10()+
    ggtitle(title)
  return(plt)
}

#plots one minus the other
plot_simple <- function(df){
  plt <- ggplot(df)+
  geom_smooth(aes(x = Nmut, y=ObservedMinusExpected), method=lm, formula = y ~ splines::bs(x, 3), colour="red")+
  #geom_smooth(aes(x = Nmut, y=ObservedUniqueMinusExpected), method='loess', formula = y ~ x, colour="orange")+
  geom_smooth(aes(x = Nmut, y=ObservedUniqueMinusExpected), method=lm, formula = y ~ splines::bs(x, 3), colour="orange")+
  geom_point(aes(x = Nmut, y=ObservedUniqueMinusExpected), alpha = 0.1, colour="orange")
  #xlim(0,250)
  #scale_x_log10()
  #scale_y_log10()
  return(plt)
}

#TEMP MOVE SOMEWHERE NEW!!!
df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/temp.tsv', sep = '\t', header=TRUE)

b <- c(1,5,10,20,30,40,50,75,100,150,1000)
ggplot(df, aes(x=Nmut_y))+
  #Draw the bars from highest to lowest --the minuses always lag one behind what is actually there
  #HIGHEST BAR--nmut strongly related
  stat_summary_bin(aes(y=nMutOncogenic), geom = 'bar', stat='mean', fill='purple', breaks = b)+
  #SECOND HIGHEST--nmut strongly related multiple
  stat_summary_bin(aes(y=nMutOncogenic - nStronglyRelatedSingleAlterations), geom = 'bar', stat='mean', fill='#31004a', breaks = b)+
  #Third highest--nmut slight related
  stat_summary_bin(aes(y=nMutOncogenic - nStronglyRelatedSingleAlterations - nStronglyRelatedMultipleMutations), geom = 'bar',  fill='#9acd32', breaks = b)+
  #fourth higest--nmut sligtly related multiple
  stat_summary_bin(aes(y=nMutOncogenic - nStronglyRelatedSingleAlterations - nStronglyRelatedMultipleMutations - nWeaklyRelatedSingleAlterations), geom = 'bar', fill='#306844', breaks = b)+
  #final bottom bar
  stat_summary_bin(aes(y=nMutOncogenic - nStronglyRelatedSingleAlterations - nStronglyRelatedMultipleMutations - nWeaklyRelatedSingleAlterations - nWeaklyRelatedMultipleMutations), geom = 'bar', fill='orange', breaks = b)+
  scale_x_log10(breaks=b)+
  stat_summary_bin(aes(y=Expected.Oncogenic), stat='mean', colour='white', breaks = b)+
  xlab('Nmut total in case')+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=12))+
  ggtitle('Genes Altered by type')
  
ggplot(head(df,6))+
  geom_bar(aes(fill=Tumor_Sample_Barcode, x=1))+
  scale_fill_manual(values = c('purple', '#31004a', '#9acd32', '#306844', 'orange', 'white'))


df$nMutWeaklyRelated
  

df <- read.table('~/Desktop/mnt/ifs/work/taylorlab/friedman/myAdjustedDataFiles/observedVsExpectedData.tsv', sep = '\t', header=TRUE)

p <-plot_simple_model_observed_vs_expected(df, 'title')

p <- plot_observed_expected_curves(df, 'Observed vs Expected (SNP+INDEL model)')
p1 <- plot_observed_expected_curves(df[df$caseMsiClass == 'Stable',], 'MSI stable')
p2 <- plot_observed_expected_curves(df[df$caseMsiClass == 'Instable',], 'MSI Instable')



################PLOTING FOR GENES
######
######
#####
#####

plot_gene_observed_expected_line_plot <- function(df, cancerType){
  df = df[df$cancer_type == cancerType,]
  df$geneMutType <- factor(df$geneMutType, levels = c('strongly_recurrent', 'weakly_recurrent', 'not_recurrent'))
  ggplot(df, aes(x = geneMutType, y=geneObsExpectedRatio))+
    geom_boxplot(fatten = NULL)+
    stat_summary(fun.y = median, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1, linetype = "solid")+
    geom_text_repel(aes(label=geneDisplayName))+
    ylim(-0.25,1)+
    ylab('Enrichment of Observed over Expected per case')+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
    ggtitle(cancerType)+
    labs(caption = 'oncogenic_mut_prob_simulation_by_gene.ipynb, compare_observed_vs_expected.R')
}

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/observedVsExpectedByGeneHypermutators.tsv', sep = '\t', header=TRUE)
Endo <- plot_gene_observed_expected_line_plot(df, 'Endometrial Cancer')
Colo <- plot_gene_observed_expected_line_plot(df, 'Colorectal Cancer')
Glio <- plot_gene_observed_expected_line_plot(df, 'Glioma')
combinedPlot <- plot_grid(Endo, Colo, Glio, ncol = 3)

ggsave('~/Desktop/plt.pdf', plot=combinedPlot,  width = 14, height = 7, units = c("in"))

#TODO TOMORROW MAKE THIS PRETTIER AND SAVE IT


df$geneDisplayName
max(df$nExpectedMuts)
















plot_observed_and_exprected_oncogenic_mut_burden <- function(df, title){
  plt <- ggplot(df[df], aes(x=reorder(Tumor_Sample_Barcode, Nmut), y=value, color=variable))+
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

plot_observed_expected_and_ratio_data <- function(df){
  thresh1 <- dim(df[df$Nmut < 25,])[1]
  thresh2 <- dim(df[df$Nmut < 100,])[1]
  plt1 <- ggplot(df, aes(x=reorder(Tumor_Sample_Barcode, Nmut)))+
    geom_point(aes(y=observedMutPerCase), colour="Red")+
    geom_point(aes(y=expectedMutPerCase), colour='Blue')+
    scale_y_log10()+
    theme(axis.title.y=element_text(size=10))+
    ylab('N Oncogenic Mutations (blue--expected, red--observed)')+
    geom_vline(aes(xintercept=thresh1))+
    geom_vline(aes(xintercept=thresh2))+
    theme(axis.title.x=element_blank(), 
          legend.title=element_blank(),
          #panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(),
          axis.line = element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank())
  
  plt2 <- ggplot(df, aes(x=reorder(Tumor_Sample_Barcode, Nmut)))+
    #geom_smooth(aes(y=ratio, group=1), method="lm")+
    geom_point(aes(y=ratio), colour='Orange')+
    ylab('Log2 ratio N expected oncogenic mutations over observed')+
    geom_hline(yintercept=0)+
    theme(axis.text.x = element_blank())+
    theme(axis.ticks.x  = element_blank())+
    geom_vline(aes(xintercept=thresh1))+
    geom_vline(aes(xintercept=thresh2))+
    geom_text(aes(x=thresh1, y=3, label='High mut burden > 25snps'), size =3)+
    geom_text(aes(x=thresh2, y=3, label='Ultra-High mut burden > 100snps'), size =3)+
    xlab('Cases sorted by Nmut')
  arrPlt <- plot_grid(plt1, plt2, align='v', rel_heights = c(.5, 1), nrow=2)
  return(arrPlt)
}

plot_observed_vs_expected_gene_level <- function(df, title, yMax=1, xMax=.3){
  plt <- ggplot(df, aes(x=expectedMutPerCase, y=observedMutPerCase, label=Hugo_Symbol))+
    geom_point(aes(size=normalCohortFrac, fill=hyperHotspotFrac), alpha = 1, pch=21)+
    geom_text(colour='#D21F3C', size=5)+
    scale_fill_viridis(option='viridis')+
    scale_size_continuous(range = c(0, 50))+
    geom_segment(aes(x=0 ,xend=xMax, y=0, yend=xMax))+
    xlim(0,xMax)+
    ylim(0,yMax)+
    ylab('N Observed Oncogenic SNPs per case ')+
    xlab('N Expected Oncogenic SNPs per case')+
    ggtitle(title)
  return(plt)
}

plot_observed_vs_expected_gene_level_simplified <- function(df, title, yMax=1, xMax=.3){
  plt <- ggplot(df, aes(x=expectedMutPerCase, y=observedMutPerCase, label=DisplayName))+
    geom_point(aes(colour=isRecurrent), alpha = 1)+
    geom_text(aes(colour=isRecurrent), size=5)+
    geom_segment(aes(x=0 ,xend=xMax, y=0, yend=xMax))+
    xlim(0,xMax)+
    ylim(0,yMax)+
    ylab('N Observed Oncogenic SNPs per case ')+
    xlab('N Expected Oncogenic SNPs per case')+
    ggtitle(title)
  plt <- plt +guides(fill=guide_legend(title="Mutated >10% in normal"))
  plt <- plt +labs(caption = "data:analyze_mutation_susceptibility_of_genes.py; plotting: compare_observed_vs_expected_oncogenic_mut_burden.R")
  return(plt)
}

plot_observed_expected_box_plot <- function(df){
  plt <- ggplot(df, aes(x=percentile, y=NOncogenicSnpsObserved))+
    geom_boxplot(fatten = NULL)+ #TODO FIX THIS TO PLOT THE MIDDLE LINE AS THE MEAN!!!
    stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1, linetype = "solid", alpha=1/2)+
    ylab('N oncogenic mutations')+
    xlab('percentile of mut burden')+
    coord_cartesian(ylim=c(0, 30))+
    ggtitle('Observed Number of Oncogenic Mutations by Percentile of Mutation Burden')
  plt<- plt+ geom_jitter(shape=16, position=position_jitter(0.1))
  
  plt <- plt + #geom_boxplot(aes(x=percentile, y=nOncogenicExpected), fatten = NULL)+ #TODO FIX THIS TO PLOT THE MIDDLE LINE AS THE MEAN!!!
         stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y.., x=percentile, y=nOncogenicSnpsExpected),
                 width = 0.75, size = 1, linetype = "solid", colour='red')
  plt <- plt + #geom_boxplot(aes(x=percentile, y=nOncogenicExpected), fatten = NULL)+ #TODO FIX THIS TO PLOT THE MIDDLE LINE AS THE MEAN!!!
    stat_summary(fun.y = median, geom = "errorbar", aes(ymax = ..y.., ymin = ..y.., x=percentile, y=nOncogenicSnpsExpected),
                 width = 0.75, size = 1, linetype = "solid", colour='orange')
  
  return(plt)
}

plot_observed_vs_expected_histogram <- function(df){
  plt <- ggplot(df)+
    stat_summary(fun.y = median, geom = "bar", aes(x=percentile2, y=NOncogenicSnpsObserved - nOncogenicSnpsExpected),
                 width = 0.75, size = 1, linetype = "solid", colour='red')
}

df <- read.table('~/Desktop/WORK/dataForLocalPlotting/simulationData.tsv', sep = '\t', header=TRUE)
p  <- plot_observed_expected_box_plot(df)
#p  <- plot_observed_vs_expected_histogram(df)
ggsave('~/Desktop/plt.pdf', p,  width = 20, height = 10, units = c("in"))

#df <- df[df$Nmut > 5,]

plt <- plot_observed_expected_and_ratio_data(df)
ggsave('~/Desktop/plt.pdf', p,  width = 20, height = 10, units = c("in"))

################GENE PLOT
df <- read.table('~/Desktop/WORK/dataForLocalPlotting/geneSpecificSimulationData_Colorectal_Cancer.tsv', sep = '\t', header=TRUE)
pltColorectal <- plot_observed_vs_expected_gene_level_simplified(df, yMax=1.7, xMax=.5, 'Observed vs Expected Tumor Suppressor Oncogenic SNP Frequencies--Colorectal')
ggsave('~/Desktop/plt.pdf', pltColorectal,  width = 15, height = 15, units = c("in"))


df <- read.table('~/Desktop/WORK/dataForLocalPlotting/geneSpecificSimulationData_Endometrial Cancer.tsv', sep = '\t', header=TRUE)
pltEndometrial <- plot_observed_vs_expected_gene_level_simplified(df, yMax=1, xMax=.5, 'Observed vs Expected Tumor Suppressor Oncogenic SNP Frequencies--Endometrial')
ggsave('~/Desktop/plt.pdf', pltEndometrial,  width = 9, height = 9, units = c("in"))

df <- read.table('~/Desktop/WORK/dataForLocalPlotting/geneSpecificSimulationData_glioma.tsv', sep = '\t', header=TRUE)
pltGlioma <- plot_observed_vs_expected_gene_level_simplified(df, yMax=.75, xMax=.3, 'Observed vs Expected Oncogenic SNP Frequencies Across 56 Hypermutated (>50Mut_mb) gliomas')
ggsave('~/Desktop/plt.pdf', pltGlioma,  width = 10, height = 10, units = c("in"))

df <- read.table('~/Desktop/WORK/dataForLocalPlotting/geneSpecificSimulationData_bladder.tsv', sep = '\t', header=TRUE)
pltBladder <- plot_observed_vs_expected_gene_level(df, yMax=1, xMax=.2, 'Observed vs Expected Oncogenic SNP Frequencies Across 9y Hypermutated (>50Mut_mb) bladders')
ggsave('~/Desktop/plt.pdf', pltBladder,  width = 15, height = 15, units = c("in"))

df <- read.table('~/Desktop/WORK/dataForLocalPlotting/geneSpecificSimulationData_melanoma.tsv', sep = '\t', header=TRUE)
pltBladder <- plot_observed_vs_expected_gene_level(df, yMax=.3, xMax=.075, 'Observed vs Expected Oncogenic SNP Frequencies Across 295 Hypermutated (>50Mut_mb) melanomas')
ggsave('~/Desktop/plt.pdf', pltBladder,  width = 15, height = 15, units = c("in"))


#################
#quick for lab pres

df <- read.table('~/Desktop/WORK/dataForLocalPlotting/indelAndSusceptibility.tsv', sep = '\t', header=TRUE)
p1<- ggplot(df, aes(x=reorder(sigName, ratio), y=ratio))+
  geom_bar(stat='identity')+
  ylab('oncogenic mutations per mutation rate')+
  ggtitle('hypermutated cancers (nmut_Mb > 50)')+
  theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(), axis.title.x=element_blank(), axis.line = element_blank())

p2 <- ggplot(df, aes(x=reorder(sigName, ratio), y=indelFrac))+
  geom_bar(stat='identity')+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
  xlab('signature class')

x <- plot_grid(
  p1, p2, nrow = 2, ncol= 1, rel_heights = c(1,.5)
)

ggsave('~/Desktop/plt.pdf', x,  width = 10, height = 10, units = c("in"))


#
#####
########
###############
########################
################################
#######################
################
##########
##

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/observedVsExpectedByGeneHypermutators.tsv', sep = '\t', header=TRUE)

ggplot(df[df$cancer_type == 'Endometrial Cancer',], aes(x = nExpectedMuts, y = nObservedMuts))+
         geom_point()+
         geom_text(aes(label=geneDisplayName))+
         geom_abline(intercept = 0)




df

