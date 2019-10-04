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



####
#CODE for overall mutation fractions
df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/observedOncFraction.tsv', sep='\t', header=TRUE)
ggplot(df[df$Nmut > 10,], aes(x = Nmut, y=fracOnc))+
  #geom_bin2d(bins = 20, aes(fill=log(..count..)))+
  geom_smooth(method='lm', formula = y ~ poly(x, 2), colour="orange")+
  geom_segment(aes(x=0, xend=700, y=.07, yend=0.07), colour='green')+
  geom_segment(aes(x=0, xend=700, y=.17, yend=0.17), colour='blue')+
  #scale_fill_viridis(direction=-1, option='plasma')+
  #geom_smooth(method=lm, formula = y ~ x)+
  #geom_point(alpha=0.1, color="blue")+
  ylim(0,1)+
  scale_x_log10()+
  ggtitle('Observed Fraction of Impact Muts That are Oncogenic')

ggplot(df, aes(x = Nmut, y=RatioGenesOncMutToGenesMut))+
  geom_smooth(method='lm', formula = y ~ poly(x, 2), colour="orange")+
  #geom_bin2d(bins = 50, aes(fill=log(..count..)))+
  #scale_fill_viridis(direction=-1, option='magma')+
  scale_x_log10()+
  ylim(0,1)+
  ggtitle('Fraction of Genes that are Oncogenically Mutated')



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
#&&&&

emptyTheme <- theme(axis.line = element_blank(),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())

#df <- read.table('~/Desktop/WORK/dataForLocalPlotting/oncogenicSusceptibilityBySignature.tsv', sep = '\t', header=TRUE)
df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/hypotheticalSignatureMutability.tsv', sep = '\t', header=TRUE)

#FIX the color columns so R can read it in appropriately
df$color = lapply(df$color, function(x) if (x != 'gray'){
  paste('#', x, sep='')}
  else{
    'gray'
  }
)

plt1<- ggplot(df, aes(x=reorder(Signature_Name, ExpectedFracOfMutsOncogenic), y=ExpectedFracOfMutsOncogenic, fill=color))+
  geom_bar(stat='identity')+
  scale_fill_identity()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  ylab('Change of a single mutation in IMPACT panel being oncogenic')+
  xlab('Signature')+
  ylim(0, 0.2)+
  ggtitle('Signature oncogenicity for SNPs')+
  emptyTheme
  #labs(caption = 'analyze_hypothetical_and_observed_mutability_by_signature.ipynb, plot_gene_mut_susceptibility.R')

plt2<- ggplot(df, aes(x=reorder(Signature_Name, oncogenicityIncludingIndels), y=oncogenicityIncludingIndels, fill=color))+
  geom_bar(stat='identity')+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  ylab('Change of a single mutation in IMPACT panel being oncogenic')+
  xlab('Signature')+
  ylim(0, 0.2)+
  ggtitle('Signature oncogenicity corrected for INDEL rates')+
  emptyTheme
  #labs(caption = 'analyze_hypothetical_and_observed_mutability_by_signature.ipynb, plot_gene_mut_susceptibility.R')

plt <- plot_grid(plt1, plt2, ncol = 2)

ggsave('~/Desktop/plot.pdf', plot=plt,  width = 15, height = 8, units = c("in"), limitsize = FALSE)

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



