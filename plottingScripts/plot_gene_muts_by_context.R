#written by Noah Friedman

#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(reshape)
library(scales)
library(forcats)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

my.rename <- function(df, old.name, new.name){ #R SUCKS!!!!!! heres my renaming function cause god forbid there would be an easy or intiuitive way to do this with R
  names(df)[names(df) == old.name] <- new.name 
  return(df)
}


plot_alteration_freqs_by_class <- function(df, title){
  p1 <- ggplot(df, aes(x = isRelatedQuadnuc, y=isRelated))+
    stat_summary(alpha=0.5)+
    ylim(0,1)+
    ylab('ccf')+
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())+
    ggtitle(title)
  p2 <- ggplot(df, aes(x = isRelatedQuadnuc))+
    geom_bar(stat='count')+
    xlab('Sig-Related QuadNuc')+
    ylab('Nmut')
  alignedPlot <- plot_grid(p1, p2, rel_heights = c(1,.5), nrow=2)
  return(alignedPlot)
}

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/mutationMotfisComp.tsv',sep = '\t', header=TRUE)
endoHot <- plot_alteration_freqs_by_class(df[df$is.a.hotspot == 'Y' & df$cancer_type == 'Endometrial Cancer',], 'Endometrial Hotspots')
endoOnc <- plot_alteration_freqs_by_class(df[df$isOncogenic == 'True' & df$cancer_type == 'Endometrial Cancer',], 'Endometrial Oncogenic')

coloHot <- plot_alteration_freqs_by_class(df[df$is.a.hotspot == 'Y' & df$cancer_type == 'Colorectal Cancer',], 'Colorectal Hotspots')
coloOnc <- plot_alteration_freqs_by_class(df[df$isOncogenic == 'True' & df$cancer_type == 'Colorectal Cancer',], 'Colorectal Oncogenic')

gliomaHot <- plot_alteration_freqs_by_class(df[df$is.a.hotspot == 'Y' & df$cancer_type == 'Glioma',], 'Glioma Hotspots')
gliomaOnc <- plot_alteration_freqs_by_class(df[df$isOncogenic == 'True' & df$cancer_type == 'Glioma',], 'Glioma Oncogenic')

alignedPlots <- plot_grid(endoOnc, endoHot, coloOnc, coloHot, gliomaOnc, gliomaHot, nrow=3, ncol=2)
title <- ggplot() + ggtitle('Mutation motif vs genes affected')
caption <- ggplot() + labs(caption = 'compare_hotspot_mutations_by_motif&plot_gene_muts_by_context.R')
plt <- plot_grid(title, alignedPlots, caption, nrow = 3, rel_heights = c(.05,1,.05))
ggsave('~/Desktop/plot.pdf', plt,  width = 6, height = 10, units = c("in"))

###NOW do this with ccf #THATS DONE BY ADJUSTING THE FUCNTION ARGUMENTS


#TEMP
df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/mutationMotfisComp.tsv',sep = '\t', header=TRUE)

p <- ggplot(df, aes(x = reorder(gene, mean), y=mean))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=7))

ggsave('~/Desktop/plot.pdf', p,  width = 10, height = 5, units = c("in"))

df$Hugo_Symbol














#ALERT CURRENTLY BEING LAZY AND ASSUMING THE Y AXIS AND FILL PARAMS are STATIC
plotFreqMutatedGenes <- function(df, xAxisDisplayParam, xAxisOrderingParam, title, plotColors, showAxisInfo=FALSE, maxColor="red"){
  df <- my.rename(df, xAxisDisplayParam, "xAxisDisplay")
  df <- my.rename(df, xAxisOrderingParam, "xAxisOrdering")
  
  df <- subset(df, !is.na(df$xAxisOrdering))
  
  plt<- ggplot(data=subset(df, !is.na(xAxisDisplay)), aes(x= reorder(xAxisDisplay, -xAxisOrdering), y=reorder(Tumor_Sample_Barcode, caseRank), fill=motif))+
    #plt<- ggplot(df, aes(x= reorder(xAxisDisplay, -xAxisOrdering), y=Tumor_Sample_Barcode, fill=fracMutated))+
    geom_tile(na.rm = TRUE)+
    #THIS CODE FILLS THE TILE
    #scale_fill_gradient(low="white",high=maxColor, limits = c(0,1))+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size=5))+
    theme(legend.position="none", axis.line = element_blank(), axis.ticks=element_blank())+
    scale_fill_manual(values=plotColors)+
    ggtitle(title)
  if(!showAxisInfo){
    plt <- plt + theme(axis.text.y=element_blank(), axis.title.y=element_blank())
  }
  else{
    plt <- plt + ylab('Tumor Sample Barcode')
  }
  if(TRUE){
    plt <- plt + theme(axis.title.x=element_blank())
  }
  return(plt)
}























#OLD CODE

plotInfreqMutatedGenes <- function(df, xAxisDisplayParam, textDisplayParam, plotColors, title, maxColor="red"){
  #df <- my.rename(df, xAxisDisplayParam, "xAxisDisplay")
  df <- my.rename(df, textDisplayParam, "textDisplay")
 
  plt <-ggplot(df,
               aes(x= uncommonOrder, y=reorder(Tumor_Sample_Barcode, caseRank), fill=motif))+
    geom_tile()+
    geom_text(aes(label=textDisplay))+
    scale_fill_manual(values=plotColors)+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size=5))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank())+
    ggtitle(title)
  if(TRUE){
    plt <- plt + theme(axis.text.y=element_blank(), axis.title.y=element_blank())
  }
  if(TRUE){
    plt <- plt + theme(axis.title.x=element_blank())
  }
}

oncoPlot_By_Ordering_Implied_By_Motif <- function(df, colors, title, sepWidth=-.055, relWidthLeft = .2){
  
  freqMutatedGenes <- plotFreqMutatedGenes(df, 'topGeneDisplayName', 'geneMutOrdering', '', colors, maxColor = "blue", showAxisInfo=TRUE)
  infreqMutatedGenes <- plotInfreqMutatedGenes(df, 'uncommonOrder', 'infrequentGeneDisplayName', colors, title,  maxColor = "blue")

  #p <- plot_grid(notTmzMotifTopGenes, notTmzMotifInfrequentGenes, tmzMotifTopGenes, tmzMotifInfrequentGenes, align='hv', nrow=1, ncol=4) 
  p <- plot_grid(freqMutatedGenes, NULL, infreqMutatedGenes,  align='hv', nrow=1, ncol=3, rel_widths = c(relWidthLeft,sepWidth,1))
  
  #return(p)
  return(p)
}



#TMZ
tmzSigGliomas <- read.table('~/Desktop/dataForLocalPlotting/tmzCasesByMotif.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
tmzSigGliomas$geneMutOrdering = sapply(tmzSigGliomas$Gene, test_order) #Placate R
tmzColors = c("#ADD8E6", "#FF0000", "Pink", "#800000")
plt <- oncoPlot_By_Ordering_Implied_By_Motif(tmzSigGliomas, tmzColors, "Oncogenic mutations in Gliomas with TMZ")
ggsave('~/Desktop/noahDogTe.pdf', plt,  width = 30, height = 8, units = c("in"))

make_bar_chart(tmzSigGliomas, 'Distribution of motifs in TMZ Gliomas', tmzColors)


#ENDOMETRIAL POLE
endometrialPoleData <- read.table('~/Desktop/dataForLocalPlotting/poleEndometrialCasesByMotif.tsv', sep = '\t', header=TRUE)

#HOTSPOT ONLY PLOT
endometrialPoleDataHotspotOnly <- endometrialPoleData[endometrialPoleData$Hotspot == 'Y', ]
endometrialPoleDataHotspotOnly$geneMutOrdering = sapply(endometrialPoleDataHotspotOnly$Gene, test_order_endometrial) #Placate R
plt <- oncoPlot_By_Ordering_Implied_By_Motif(endometrialPoleDataHotspotOnly, poleColors, "Oncogenic mutations in POLE Endometrials by Mutation Motif",sepWidth=-.045, relWidthLeft = .175)
ggsave('~/Desktop/noahDogTeo.pdf', plt,  width = 40, height = 8, units = c("in"))


endometrialPoleData$geneMutOrdering = sapply(endometrialPoleData$Gene, test_order_endometrial) #Placate R
poleColors = c("#00FFFF", "gray", "#017455", "#ADFF2F", "#00DD5D") #ORDER: agingonly, other, poleAndAging, PoleOPnly, Sig14Only
plt <- oncoPlot_By_Ordering_Implied_By_Motif(endometrialPoleData, poleColors, "Oncogenic mutations in POLE Endometrials by Mutation Motif",sepWidth=-.045, relWidthLeft = .175)
ggsave('~/Desktop/noahDogTeo.pdf', plt,  width = 40, height = 8, units = c("in"))

endometrialPoleAllMutsData <- read.table('~/Desktop/dataForLocalPlotting/poleEndometrialAllMutsForBarchart.tsv', sep = '\t', header=TRUE)
#BARCHARTS For Context
make_bar_chart(endometrialPoleAllMutsData, 'Distribution of motifs of all mutations in POLE hypermutators', poleColors)

endometrialPoleDataTop5 = endometrialPoleData[endometrialPoleData$Gene %in% c("TP53", "PTEN", 'PIK3CA', 'ARID1A', 'POLE'),]

make_bar_chart(endometrialPoleDataTop5[endometrialPoleDataTop5$Hotspot == 'Y', ], 'Distribution of motifs of top mutated genes POLE HOTSPOTS', poleColors)



