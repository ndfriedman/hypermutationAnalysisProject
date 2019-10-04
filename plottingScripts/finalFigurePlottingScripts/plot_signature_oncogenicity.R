#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)

#TODO MAKE THIS CODE BETTER I WAS RUSHED WHEN I MADE IT

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
  plt<- #ggplot(df, aes(x=Nmut_Mb, y=oncogenicMutsPerCase, colour=dominantSigAdj)) +
    ggplot(df, aes(x=Nmut_Mb, y=oncogenicMutsPerCase, group =dominantSigAdj, colour=sigDisplayLabel)) +
    geom_smooth(method='loess', formula = y ~ x, alpha=0.25)+
    coord_cartesian(xlim = c(10,1000))+
    scale_x_log10()+
    ggtitle('Observed Oncogenicty\n of Signatures')+
    ylab('N Observed Oncogenic Mutations')+
    scale_colour_manual(values = palette, name = 'Dominant\nSignature')+
    emptyTheme
    return(plt)
}


emptyTheme <- theme(axis.line = element_blank(),
                    axis.ticks.x = element_blank(),
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

pltOncSnp<- ggplot(df, aes(x=reorder(Signature_Name, ExpectedFracOfMutsOncogenic), y=ExpectedFracOfMutsOncogenic, fill=color))+
  geom_bar(stat='identity')+
  scale_fill_identity()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  ylab('Chance of oncogenic SNP')+
  xlab('Signature')+
  ylim(0, 0.2)+
  ggtitle('Signature oncogenicity for SNPs')+
  emptyTheme

pltOncAll <- ggplot(df, aes(x=reorder(Signature_Name, oncogenicityIncludingIndels), y=oncogenicityIncludingIndels, fill=color))+
  geom_bar(stat='identity')+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  ylab('Chance of oncogenic mutation')+
  xlab('Signature')+
  ylim(0, 0.2)+
  ggtitle('Signature oncogenicity corrected for INDEL rates')+
  emptyTheme

df2 <- read.table('~/Desktop/WORK/dataForLocalPlotting/observedOncogenicDataBySig.tsv', sep='\t', header=TRUE)
pltObservedOncogenicity <- plot_observed_mutations_by_signature_aetiology(df2)

alignedOncogenicities <- plot_grid(pltOncSnp, pltOncAll, nrow=2)
fullPlot <- plot_grid(alignedOncogenicities, pltObservedOncogenicity, ncol=2)

ggsave('~/Desktop/plot.pdf', plot=fullPlot,  width = 9, height = 6, units = c("in"))










