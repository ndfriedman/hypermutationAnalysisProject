#written by noah friedman
#a template for R scripts for plotting
#A script that plots the pan cohort numbers of activating mutations by mutation burden and the ratio of activating mutations to total number of mutations

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)

graph_mutation_ratio_prevalences <- function(df, xValParam, yColParam1, yColParam2, xlabel, yMax1, yMax2, nbins=50){
  
  df <- my.rename(df, xValParam, "xVal")
  df <- my.rename(df, yColParam1, "yCol1") 
  df <- my.rename(df, yColParam2, "yCol2")
  plt<-ggplot()+
    stat_summary_bin(data=df, aes(x=xVal, y=yCol1), fun.y = "mean", geom = "bar",bins = nbins)+
    scale_x_log10()+
    coord_cartesian(ylim=c(0, yMax1))+
    xlab(xlabel)+
    ylab(yColParam1)+
    stat_summary_bin(data=df, aes(x=xVal, y = yCol2*(yMax1/yMax2)), fun.y = "mean", geom = "path",bins = nbins, colour="red")
  plt<-plt+scale_y_continuous(sec.axis = sec_axis(~ . *(yMax2/yMax1) , name = "N activating mut/Nmut Ratio"))
  plt<-plt+ggtitle("Total mutation burden vs number of activating mutations")
  return(plt)
}


graph_mutation_ratio_prevalences2 <- function(df, xValParam, yColParam1, xlabel, yMax1, yMax2, title, nbins=50){
  
  df <- my.rename(df, xValParam, "xVal")
  df <- my.rename(df, yColParam1, "yCol1") 
  plt<-ggplot()+
    stat_summary_bin(data=df, aes(x=xVal, y=yCol1), fun.y = "mean", geom = "bar",bins = nbins)+
    scale_x_log10()+
    coord_cartesian(ylim=c(0, yMax1))+
    xlab(xlabel)+
    ylab(yColParam1)
    plt<-plt+ggtitle(title)
  return(plt)
}


if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")


plot_mut_rate_by_mutation_type <- function(df, title=''){
  plt <- ggplot(data=df, aes(x=Nmut_Mb_Case))+
    geom_smooth(aes(y=Nmut_Related, colour='related'), method='loess', formula = y ~ x)+
    geom_smooth(aes(y=Nmut_Non.Related, colour='unrelated'), method='loess', formula = y ~ x)+
    #scale_x_log10()+
    #scale_y_log10(breaks=c(1,5,10,25,100))+
    geom_point(aes(y=Nmut_Related, colour='related'), alpha = 0.05)+
    geom_point(aes(y=Nmut_Non.Related, colour='unrelated'), alpha = 0.05)+
    scale_colour_manual(name="GeneType",
                      values=c(unrelated="orange", related="purple"))+
  
    ylab('N oncogenic mutations in case')+
    xlab('Nmut_Mb')+
    ggtitle(title)+
    labs(caption = "data:compare_hypermutator_vs_normal_gene_mut_characteristics.py; plotting: activating_mut_vs_nmut_plotting.R")
  
  return(plt)
}

#COMPARE MSI/OTHER STUFF
plot_mut_rate_by_mutation_type_msi_focus <- function(df, title=''){
  plt <- ggplot()+
    geom_smooth(data=df[df$caseMsiClass != 'Instable',], aes(x=Nmut_Mb_Case, y=Nmut_Uniq_Related, colour='related_MSI_STABLE'), method='loess', formula = y ~ x)+
    geom_smooth(data=df[df$caseMsiClass != 'Instable',], aes(x=Nmut_Mb_Case, y=Nmut_Uniq_Non.Related, colour='unrelated_MSI_STABLE'), method='loess', formula = y ~ x)+
    geom_smooth(data=df[df$caseMsiClass == 'Instable',], aes(x=Nmut_Mb_Case, y=Nmut_Uniq_Related, colour='related_MSI_INSTABLE'), method='loess', formula = y ~ x)+
    geom_smooth(data=df[df$caseMsiClass == 'Instable',], aes(x=Nmut_Mb_Case, y=Nmut_Uniq_Non.Related, colour='unrelated_MSI_INSTABLE'), method='loess', formula = y ~ x)+
    scale_x_log10()+
    geom_point(data=df[df$caseMsiClass != 'Instable',], aes(x=Nmut_Mb_Case, y=Nmut_Uniq_Related, colour='related_MSI_STABLE'), alpha=0.05)+
    geom_point(data=df[df$caseMsiClass != 'Instable',], aes(x=Nmut_Mb_Case, y=Nmut_Uniq_Non.Related, colour='unrelated_MSI_STABLE'), alpha=0.05)+
    geom_point(data=df[df$caseMsiClass == 'Instable',], aes(x=Nmut_Mb_Case, y=Nmut_Uniq_Related, colour='related_MSI_INSTABLE'), alpha=0.1)+
    geom_point(data=df[df$caseMsiClass == 'Instable',], aes(x=Nmut_Mb_Case, y=Nmut_Uniq_Non.Related, colour='unrelated_MSI_INSTABLE'), alpha=0.1)+
  
    scale_colour_manual(name="GeneType",
                        values=c(unrelated_MSI_STABLE="orange", related_MSI_STABLE="#DDA0DD", related_MSI_INSTABLE="purple", unrelated_MSI_INSTABLE="#FF4500"))+
    
    ylab('N distinct genes oncogenically mutated')+
    xlab('Nmut_Mb')+
    ggtitle(title)+
    labs(caption = "data:compare_hypermutator_vs_normal_gene_mut_characteristics.py; plotting: activating_mut_vs_nmut_plotting.R")
  
  return(plt)
}

#COMPARE TUMOR SUPRESSORS/OTHER STUFF
plot_mut_rate_by_mutation_type_exteneded <- function(df, title=''){
  plt <- ggplot(data=df, aes(x=Nmut_Mb_Case))+
    geom_smooth(aes(y=Nmut_Related_Oncogene, colour='relatedOncogene'), method='loess', formula = y ~ x)+
    geom_smooth(aes(y=Nmut_Non.Related_Oncogene, colour='unrelatedOncogene'), method='loess', formula = y ~ x)+
    geom_smooth(aes(y=Nmut_Related_Tumor_Supressor, colour='relatedTumorSupressor'), method='loess', formula = y ~ x)+
    geom_smooth(aes(y=Nmut_Non.Related_Tumor_Supressor, colour='unrelatedTumorSupressor'), method='loess', formula = y ~ x)+
    scale_x_log10()+
    scale_colour_manual(name="GeneType",
                        values=c(unrelatedOncogene="pink", relatedOncogene="red", relatedTumorSupressor="blue", unrelatedTumorSupressor="#99CCFF"))+
    
    ylab('N oncogenic mutations in case')+
    xlab('Nmut_Mb')+
    ggtitle(title)+
    labs(caption = "data:compare_hypermutator_vs_normal_gene_mut_characteristics.py; plotting: activating_mut_vs_nmut_plotting.R")
  return(plt)
}

df <- read.table('~/Desktop/WORK/dataForLocalPlotting/oncogenicMutCount.tsv',sep = '\t', header=TRUE)

df

dfSmall <- df[1:4000,]

plt <- plot_mut_rate_by_mutation_type(dfSmall, title='Mutations by type')
plt <- plot_mut_rate_by_mutation_type_exteneded(df, title='Mutations by type')

plt <- plot_mut_rate_by_mutation_type_msi_focus(df, title='Mutations by type')

ggsave('~/Desktop/testNoah.pdf', plot=plt, width = 6, height = 4, units = c("in"))


























df <- read.table('/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myAdjustedDataFiles/mafSummaryData/mafOncogenicMutSummaryData.tsv',sep = '\t', header=TRUE)

#df$hotspotMutToNmutRatio df$oncogenicMutToNmutRatio df$nHotspotMuts df$nOncogenicMuts
dfAdj <- df[df$Nmut < 500 & df$Nmut>50,]

plt <- graph_mutation_ratio_prevalences(df,
                                 xValParam="Nmut", 
                                 yColParam1="nActivatingMutations", 
                                 yColParam2="activatingMutToNmutRatio", 
                                 xlabel="Number of Mutations (log-scale)", 
                                 yMax1=max(df$nActivatingMutations), 
                                 yMax2=max(df$activatingMutToNmutRatio),
                                 nbins=50)

ggsave('~/Desktop/testNoah.pdf',plot=plt)


####################################

dfEndometrial <- read.table('~/Desktop/WORK/dataForLocalPlotting/endometrialCNAData.tsv', sep = '\t', header=TRUE)
dfBladder <- read.table('~/Desktop/WORK/dataForLocalPlotting/bladderCNAData.tsv', sep = '\t', header=TRUE)
dfColorectal <- read.table('~/Desktop/WORK/dataForLocalPlotting/colorectalCNAData.tsv', sep = '\t', header=TRUE)
dfGlioma <- read.table('~/Desktop/WORK/dataForLocalPlotting/gliomaCNAData.tsv', sep = '\t', header=TRUE)
dfMelanoma <- read.table('~/Desktop/WORK/dataForLocalPlotting/melanomaCNAData.tsv', sep = '\t', header=TRUE)


pltEndometrial <- graph_mutation_ratio_prevalences2(dfEndometrial,
                                        xValParam="Nmut_Mb", 
                                        yColParam1="FGA", 
                                        xlabel="Nmut_Mb (log-scale)", 
                                        yMax1=.25, 
                                        yMax2=max(df$Nmut),
                                        title="FGA by mutation burden: Endometrial",
                                        nbins=10)

pltBladder <- graph_mutation_ratio_prevalences2(dfBladder,
                                                    xValParam="Nmut_Mb", 
                                                    yColParam1="FGA", 
                                                    xlabel="Nmut_Mb (log-scale)", 
                                                    yMax1=.25, 
                                                    yMax2=max(df$Nmut),
                                                    title="FGA by mutation burden: Bladder",
                                                    nbins=10)

pltColorectal <- graph_mutation_ratio_prevalences2(dfColorectal,
                                                    xValParam="Nmut_Mb", 
                                                    yColParam1="FGA", 
                                                    xlabel="Nmut_Mb (log-scale)", 
                                                    yMax1=.25, 
                                                    yMax2=max(df$Nmut),
                                                    title="FGA by mutation burden: Colorectal",
                                                    nbins=10)

pltGlioma <- graph_mutation_ratio_prevalences2(dfGlioma,
                                                    xValParam="Nmut_Mb", 
                                                    yColParam1="FGA", 
                                                    xlabel="Nmut_Mb (log-scale)", 
                                                    yMax1=.25, 
                                                    yMax2=max(df$Nmut),
                                                    title="FGA by mutation burden: Glioma",
                                                    nbins=10)

pltMelanoma <- graph_mutation_ratio_prevalences2(dfMelanoma,
                                                    xValParam="Nmut_Mb", 
                                                    yColParam1="FGA", 
                                                    xlabel="Nmut_Mb (log-scale)", 
                                                    yMax1=.25, 
                                                    yMax2=max(df$Nmut),
                                                    title="FGA by mutation burden: Melanoma",
                                                    nbins=10)


plt2 <-plot_grid(
  pltBladder, pltGlioma, pltMelanoma, pltEndometrial, pltColorectal, 
  align='hv', nrow=2, ncol=3
)

title <- ggdraw() + draw_label("Fraction Genome Altered (FGA) by mutation burden", fontface='bold')
finalPlot <- plot_grid(title, plt2, ncol=1, rel_heights=c(0.1, 1))

ggsave('~/Desktop/testNoah.pdf', plot=finalPlot, width = 20, height = 15, units = c("in"))

#
###
########
############
###############
##################
########################
#############################
#################################
#############################
#########################
######################
#################
#############
#########
#####
#

emptyTheme <- theme(axis.line = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/relatedGeneFracByCancerType.tsv', sep = '\t', header=TRUE)
p <- ggplot(df, aes(x = Nmut_Mb, y = relatedFrac, group=cancer_type, colour=cancer_type))+
  geom_smooth(formula = y ~ x, method='lm', se=FALSE, alpha=0.75)+
  ylim(0,1)+
  ylab('Fraction of mutations in Related Genes')+
  emptyTheme+
  scale_colour_brewer(palette = "Set1")+
  ggtitle('Mutations in Cancer Type Related Genes')+
  labs(caption='activating_mut_vs_nmut_plotting.R, hotspot_and_motif_analysis.ipynb')

ggsave('~/Desktop/plot.pdf', plot=p, width = 5, height = 4, units = c("in"))





#
####
#########
############
#################
#####################
###############
###########
#########
######
##
#

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/geneMutTypeSummary.tsv', sep = '\t', header=TRUE)

plt <- ggplot(df, aes(x=tmb))+
  geom_smooth(aes(y=nRelated, colour='N mutations in DNDS\n strongly significant genes'))+
  geom_smooth(aes(y=nHyperMuts, colour='N mutations in DNDS\n hypermutation only significant genes'))+
  geom_smooth(aes(y=nWeakMuts, colour='N mutations in DNDS\n weakly significant genes'))+
  geom_smooth(aes(y=nOncogenic - nRelated - nHyperMuts - nWeakMuts, colour='N mutations in\nDNDS not-significant genes'))+
  scale_x_log10()+
  ylab('N mutations')+
  xlab('Log TMB')+
  ggtitle('Oncogenic mutations by gene mutated and TMB')


ggsave('~/Desktop/plot.pdf', plt,  width = 3, height = 4, units = c("in"))


emptyTheme <- theme(axis.line = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(), 
                    axis.ticks = element_blank(),
                    axis.text= element_blank())

