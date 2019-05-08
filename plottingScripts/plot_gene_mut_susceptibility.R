#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)
library(ggrepel)
library(RColorBrewer) 

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

layer_gene_plot <- function(df){
  plt1 <- ggplot(df, aes(x=reorder(gene, -orderingVal), y=nDistinctOncogenicSnps))+
    geom_bar(stat='identity')+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))+
    ggtitle('Number of distinct oncogenic mutations observed in genes')
  plt2 <- ggplot(df, aes(x=reorder(gene, -orderingVal), y=cdsLength))+
    geom_bar(stat='identity')+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))+
    labs(caption = "data:analyze_mutation_susceptibility.py; plotting: plot_gene_mut_susceptibility.R")
  alignedPlot <- plot_grid(
    plt1, plt2,
    align='hv', nrow=2, ncol=1
  )
  return(alignedPlot)
}

quadNucMutSusceptibilityPlot <- function(df, title=''){
  plt <- ggplot(df, aes(x=reorder(QuadNuc, orderingVal), y=frac_oncogenic, fill=factor(orderingVal)))+
    geom_bar(stat = 'identity')+
    scale_fill_manual(values=c('#1EBFF0', '#050708', '#E62725', '#CBCACB', '#A1CF64', '#EDC8C5'))+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))+
    ggtitle(title)+
    xlab('Mutation and Motif')+
    theme(legend.position="none")+
    labs(caption = "data=mutation_modeling_workspace.py--from /ifs/work/taylorlab/friedman/myAdjustedDataFiles/impact341SimulatedMutationDataSummary.tsv; plot:plot_gene_mut_susceptibility.R")
  return(plt)
}

plot_expected_mutations_by_signature_aetiology<- function(df){
  plt<- ggplot(df, aes(x=Nmut_Mb, y=Nmut_Expected, colour=Signature, alpha=colorOpaque)) +
    geom_line()+
    scale_x_log10()+
    ggtitle('Oncogenicty of Signatures')+
    ylab('N Expected Oncogenic Mutations')+
    theme(legend.position="none")+
    labs(caption="data:compare_hypothetical_mutability_by_signature.py, plot_gene_mut_susceptibility.R")
  return(plt)
}

plot_observed_mutations_by_signature_aetiology<- function(df){
 palette = c( 
  "#00DFFF", #age is 
  "#FF0000", #apobec is red
  "#FF1493", #BRCA sig deep pink
  "#267574", #mmr is blue-green
  "gray", #other is gray
  "#ADFF2F", #Pole is
  "#551A8B", #Pole_MMR
  "#FFA500", #smoking is orange
  "purple", #tmz
  "#FFF600" #UV is
  )
  plt<- ggplot(df, aes(x=Nmut_Mb, y=oncogenicMutsPerCase, colour=dominantSigAdj)) +
    #ggplot(df, aes(x=Nmut_Mb, y=oncogenicMutsPerCase, group =dominantSigAdj, colour=sigDisplayLabel)) +
    geom_smooth(method='loess', formula = y ~ x, alpha=0.25)+
    coord_cartesian(xlim = c(10,1000))+
    scale_x_log10()+
    ggtitle('Oncogenicty of Signatures')+
    ylab('N Observed Oncogenic Mutations')
    #scale_colour_manual(values = palette)
    #theme(legend.position="none")+
    #labs(caption="data:compare_hypothetical_mutability_by_signature.py, plot_gene_mut_susceptibility.R")
  return(plt)
}


###CURRENT CODE
#PLOT expected curves
df <- read.table('~/Desktop/WORK/dataForLocalPlotting/expectedOncogenicDataBySig.tsv', sep = '\t', header=TRUE)
plt <- plot_expected_mutations_by_signature_aetiology(df)
ggsave('~/Desktop/plot.pdf', plot=plt,  width = 3, height = 4, units = c("in"), limitsize = FALSE)

#PLOT observed curves
df <- read.table('~/Desktop/WORK/dataForLocalPlotting/observedOncogenicDataBySig.tsv', sep='\t', header=TRUE)
plt <- plot_observed_mutations_by_signature_aetiology(df)
ggsave('~/Desktop/plot.pdf', plot=plt,  width = 5, height = 4, units = c("in"), limitsize = FALSE)








##########OLDER CODE
#####STUFF WITH FUNCTIONS


#oncogenic mut susceptibility based on quad nuc motif of a mutation
df <- read.table('~/Desktop/WORK/dataForLocalPlotting/quadNucOncogenicInfo.tsv', sep = '\t', header=TRUE)
plt <- quadNucMutSusceptibilityPlot(df, title='Percentage of mutations at a quad nuc that are oncogenic in impact')
ggsave('~/Desktop/plot.pdf', plot=plt,  width = 20, height = 6, units = c("in"), limitsize = FALSE)
#
#observed mut susceptibility
df <- read.table('~/Desktop/WORK/dataForLocalPlotting/quadNucOncogenicInfo_Observed.tsv', sep = '\t', header=TRUE)
plt <- quadNucMutSusceptibilityPlot(df, title='Observed fraction of mutations that are oncogenic at quadnuc')
barP <- ggplot(df, aes(x=reorder(QuadNuc, orderingVal), y=Nmut))+
  geom_bar(stat='identity')+
  get_empty_theme()+
  theme(axis.title.y=element_text(hjust = .2, size=7, face="bold"))+
  ylab('N mutations at motif in impact')
alginedPlot <- plot_grid(
  plt, barP, 
  align='v',nrow=2, ncol=1, rel_heights = c(1,.2)
)
ggsave('~/Desktop/plot.pdf', plot=alginedPlot,  width = 20, height = 6, units = c("in"), limitsize = FALSE)
  
##########SIGNATURE BASED MUT SUSCEPTIBILITY  

df <- read.table('~/Desktop/WORK/dataForLocalPlotting/oncogenicSusceptibilityBySignature.tsv', sep = '\t', header=TRUE)
plt<- ggplot(df, aes(x=reorder(Signature_Name, ExpectedFracOfMutsOncogenic), y=ExpectedFracOfMutsOncogenic))+
  geom_bar(stat='identity')+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  ylab('Change of a single mutation in IMPACT panel being oncogenic')+
  xlab('Signature')+
  ggtitle('Signature propensity to cause oncogenic mutations in IMPACT 341')

ggsave('~/Desktop/plot.pdf', plot=plt,  width = 10, height = 10, units = c("in"), limitsize = FALSE)


######################  
#MUT susceptibility as number of distinct oncogenic mutations/size of impact region
df <- read.table('~/Desktop/WORK/dataForLocalPlotting/geneOncogenicInfo.tsv',sep = '\t', header=TRUE)
plt <- layer_gene_plot(df)
ggsave('~/Desktop/plot.pdf', plot=plt,  width = 75, height = 8, units = c("in"), limitsize = FALSE)


#Truncating mutations vs Possible stop codons
df1 <- read.table('~/Desktop/WORK/dataForLocalPlotting/truncatingInfoPole.tsv',sep = '\t', header=TRUE)

dfZ <- df1[df1$NObsercedTruncAtMotif > 0,]

p <- truncating_gene_plot(dfZ, title='Possible Stop Gain Mutations Introduced by POLE vs Observed Stop Gain Mutations Endometrial POLE (n=23)')
ggsave('~/Desktop/plot.pdf', plot=p,  width = 45, height = 8, units = c("in"), limitsize = FALSE)

ggsave('~/Desktop/plot.pdf', plot=p,  width = 45, height = 45, units = c("in"), limitsize = FALSE)



