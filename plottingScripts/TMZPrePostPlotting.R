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

#MAKES an up and down bar plot of the most frequently mutated genes in TMZ gliomas, colored by mut freq catergory
barPlotGeneFreq <- function(df){
  plt <- ggplot(df)+
    #The first bar of the predominant signature in the positive direction
    geom_bar(aes(x = Tumor_Sample_Barcode, y=NPreTMZOncogenicMuts, fill=class), stat="identity")+
    #The second bar of the second predominant signature in the negative direction
    geom_bar(aes(x = Tumor_Sample_Barcode, y=-NPostTMZOncogenicMuts, fill=class), stat="identity")+
    #scale_fill_manual()
    theme(legend.text=element_text(size=3), legend.title=element_text(size=5))+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=5))
}

my.rename <- function(df, old.name, new.name){ #R SUCKS!!!!!! heres my renaming function cause god forbid there would be an easy or intiuitive way to do this with R
  names(df)[names(df) == old.name] <- new.name 
  return(df)
}

#ALERT CURRENTLY BEING LAZY AND ASSUMING THE Y AXIS AND FILL PARAMS are STATIC
plotFreqMutatedGenes <- function(df, xAxisDisplayParam, xAxisOrderingParam, title, showAxisInfo=FALSE, maxColor="red"){
  df <- my.rename(df, xAxisDisplayParam, "xAxisDisplay")
  df <- my.rename(df, xAxisOrderingParam, "xAxisOrdering")
  
  df <- subset(df, !is.na(df$xAxisOrdering))
  
  plt<- ggplot(data=subset(df, !is.na(xAxisDisplay)), aes(x= reorder(xAxisDisplay, -xAxisOrdering), y=reorder(Tumor_Sample_Barcode, caseRank), fill=fracMutated))+
    #plt<- ggplot(df, aes(x= reorder(xAxisDisplay, -xAxisOrdering), y=Tumor_Sample_Barcode, fill=fracMutated))+
    geom_tile(na.rm = TRUE)+
    #THIS CODE FILLS THE TILE
    scale_fill_gradient(low="white",high=maxColor, limits = c(0,.6))+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size=5))+
    theme(legend.position="none", axis.line = element_blank(), axis.ticks=element_blank())+
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

plotInfreqMutatedGenes <- function(df, xAxisDisplayParam, textDisplayParam, maxColor="red"){
  df <- my.rename(df, xAxisDisplayParam, "xAxisDisplay")
  df <- my.rename(df, textDisplayParam, "textDisplay")
  
  df <- subset(df, !is.na(df$xAxisDisplay))
  
  plt <-ggplot(df,
               aes(x= reorder(xAxisDisplay, -fracMutated), y=reorder(Tumor_Sample_Barcode, caseRank), fill=fracMutated))+
    geom_tile()+
    geom_text(aes(label=textDisplay))+
    #THIS CODE FILLS THE TILE
    scale_fill_gradient(low="white",high=maxColor, limits = c(0, .6))+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size=5))+
    theme(legend.position="none")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank())
  if(TRUE){
    plt <- plt + theme(axis.text.y=element_blank(), axis.title.y=element_blank())
  }
  if(TRUE){
    plt <- plt + theme(axis.title.x=element_blank())
  }
}

oncoPlotPrePostTMZ <- function(df){
  
  preTMZTopGenes <- plotFreqMutatedGenes(df, 'topGliomaGeneDisplayNamePreTMZ', 'testOrderPre', 'Mutations before TMZ', maxColor = "blue", showAxisInfo=TRUE)
  preTMZInfrequentGenes <- plotInfreqMutatedGenes(df, 'preTMZUncommonOrder', 'Gene', maxColor = "blue")
  postTMZTopGenes <- plotFreqMutatedGenes(df, 'topGliomaGeneDisplayNamePostTMZ', 'testOrderPost', 'Mutations after TMZ')
  postTMZInfrequentGenes <- plotInfreqMutatedGenes(df, 'postTMZUncommonOrder', 'Gene')
  
  #USE a bunch of NULL plots with negative relative widths to make everything align appropriately
  p <- plot_grid(
    preTMZTopGenes+ theme(plot.margin=margin(r=-2.8,unit="in")),
    NULL,
    preTMZInfrequentGenes+ theme(plot.margin=margin(r=-2.8,unit="in")),
    NULL,
    postTMZTopGenes+ theme(plot.margin = unit(c(0, 0, 0, 0), "in")),
    NULL,
    postTMZInfrequentGenes + theme(plot.margin = unit(c(0, 0, 0, 0), "in")),
    align='hv', nrow=1, ncol=7, rel_widths = c(.2, -.08, .1, -.025, .125, -.075, 1)
  )
  return(p)      
}

#TODO fix #A BRAIN DEAD COLUMN ORDERING TO PLACATE THE DIABOLICAL CAPRICES OF R
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

tmzPrePost <- read.table('~/Desktop/dataForLocalPlotting/tmzMultipleSampleGenes.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util

tmzPrePost$testOrderPre = sapply(tmzPrePost$topGliomaGeneDisplayNamePreTMZ, test_order) #Placate R
tmzPrePost$testOrderPost = sapply(tmzPrePost$topGliomaGeneDisplayNamePostTMZ, test_order)
plt <- oncoPlotPrePostTMZ(tmzPrePost)

ggsave('~/Desktop/noahDogTe.pdf', plt,  width = 30, height = 5, units = c("in"))

  


####################################
df <- read.table('~/Desktop/WORK/dataForLocalPlotting/gliomaMMRGenes.tsv', sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util


plt1 <- ggplot(df, aes(x=orderingVal, y=reorder(Tumor_Sample_Barcode, genomeDoubled)))+
  geom_text(aes(label=Hugo_Symbol))+
  geom_tile(aes(fill=Oncogenic), alpha=0.3)+
  theme(axis.text.x = element_blank())+
  theme(axis.ticks.x  = element_blank())+
  theme(axis.text.y = element_blank())+
  ylab('')+
  theme(axis.ticks.y  = element_blank())

plt2 <- ggplot(df, aes(x=1, y=reorder(Tumor_Sample_Barcode, genomeDoubled)))+
  geom_tile(aes(fill=genomeDoubled))+
  theme(legend.position="none")+
  theme(axis.ticks.x  = element_blank())+
  theme(axis.text.x = element_blank())

alignedPlt <- plot_grid(plt2, plt1, ncol = 2, rel_widths = c(.2,1)) 

ggsave('~/Desktop/noahDogTe.pdf', alignedPlt,  width = 15, height = 15, units = c("in"))





