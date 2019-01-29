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


test_order <- function(x){
  if(x == 'TERT'){
    return(10)
  }
  if(x == 'TP53'){
    return(9)
  }
  if(x == 'IDH1'){
    return(8)
  }
  if(x == 'PTEN'){
    return(7)
  }
  if(x == 'ATRX'){
    return(6)
  }
  else{
    return(NA)
  }
}

test_order_endometrial <- function(x){
  if(x == 'POLE'){
    return(11)
  }
  if(x == 'TP53'){
    return(10)
  }
  if(x == 'PIK3CA'){
    return(9)
  }
  if(x == 'PTEN'){
    return(8)
  }
  if(x == 'ARID1A'){
    return(7)
  }
  if(x == 'PIK3R1'){
    return(6)
  }
  if(x == 'KRAS'){
    return(5)
  }
  if(x == 'CTNNB1'){
    return(4)
  }
  if(x == 'CTCF'){
    return(3)
  }
  if(x == 'PPP2R1A'){
    return(2)
  }
  if(x == 'FBXW7'){
    return(1)
  }
  else{
    return(NA)
  }
}

make_bar_chart <- function(df, title, colors){
  ggplot(df, aes(x="", fill=motif))+
    geom_bar(width = 1, stat = "count")+
    coord_polar("y", start=0)+
    ggtitle(title)+
    scale_fill_manual(values=colors)+
    theme_minimal()
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



