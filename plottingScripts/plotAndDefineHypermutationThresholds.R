#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)

library("mclust")
library(mclust, quietly=TRUE)
library(gridBase)
if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

plot_distribution <- function(df, binW = 0.5, title=''){
  plt<- ggplot(df, aes(x=log(Nmut_Mb), fill=factor(hypermutantClassification)))+
    geom_histogram(aes(y=..count../sum(..count..)), binwidth=binW)+
    geom_density(alpha=.2, fill="#FF6666")+
    scale_x_continuous(breaks=c(1,2,2.5,3,3.5,4,4.5,8))+
    theme(axis.text.x=element_text(size=3, face="bold"))+
    ggtitle(title)
  return(plt)
}

fit_and_analyze_dist <- function(dataFrame){
  dataFrame = dataFrame[!is.na(dataFrame$Nmut_Mb),]
  Y <- dataFrame$Nmut_Mb
  obj = densityMclust(Y, G=2)
  
  df <- data.frame(c(obj$data), c(obj$classification), dataFrame$Tumor_Sample_Barcode)
  
  colN <- c("Nmut_Mb","Classification", 'Tumor_Sample_Barcode')
  colnames(df) <- colN
  
  dimClass1 <- dim(df[df$Classification == 1,])[1]
  dimClass2 <- dim(df[df$Classification == 2,])[1]
  smallerClassSize <- min(dimClass1, dimClass2)/(dimClass1 + dimClass2)
  smallerClassData <- df[df$Classification == 2,]$Nmut_Mb #THIS IS THE HYPERMUTATORS THE GROUP WITH A SMALL CLASS SIZE
  largerClassData <- df[df$Classification == 1,]$Nmut_Mb #THIS IS THE NON HYPERMUTATORS THE GROUP WITH A LARGE CLASS SIZE
   
  #IF THE SMALLER CLASS IS TOO BIG THAT IS AN ISSUE AND WE NEED TO SIGNAL AN ERROR
  
  indeterminateUpperBound <- max(quantile(smallerClassData, c(.2)), 20)
  indeterminateLowerBound <- min(20, quantile(smallerClassData, c(0)))
  
  df <- mutate(df,
         hypermutantClassification = ifelse(Nmut_Mb>=indeterminateUpperBound, "Hypermutated",
                                    ifelse(Nmut_Mb>=indeterminateLowerBound, "Indeterminate",
                                    ifelse(Nmut_Mb<indeterminateLowerBound, "Normal", "No_Nmut_Mb_Info")))
                                            )
  df$HypermutantThresh = indeterminateUpperBound
  df$NormalThresh = indeterminateLowerBound
  
  return(list(df, smallerClassData, largerClassData))
}


#sigsData <- read.table('~/Desktop/WORK/dataForLocalPlotting/sigsWithCType.tsv', sep='\t', header=TRUE)
sigsData <- read.table('~/Desktop/mnt/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/tmbInfo.tsv', sep='\t', header=TRUE)

l <- list()
i <- 1
cancerTypes <- unique(sigsData$cancer_type)
minNSamples <- 250
for (cancerType in cancerTypes){
  analyzeD = sigsData[sigsData$cancer_type == cancerType,]
  if(dim(analyzeD)[[1]] >= minNSamples){
    info <- fit_and_analyze_dist(analyzeD)
    plotData <- info[[1]]
    dataLarge<- info[[2]]
    dataSmall <- info[[3]]
    p <- plot_distribution(plotData, title=cancerType)
    l[[i]] <- p
    i <- i + 1
    
    if(nchar(cancerType) > 1){
      path <- '~/Desktop/mnt/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/hypermutationStatusIds/'
      cancerTypeAdj <- gsub(" ", "_", cancerType)
      fullPath <- paste(path, cancerTypeAdj, '.tsv', sep='')
      write.table(plotData, file=fullPath, sep='\t')
      }
    }
}

length(l)

p <- plot_grid(l[[1]], l[[2]], l[[3]], l[[4]], l[[5]],
          l[[6]], l[[7]], l[[8]], l[[9]], l[[10]],
          l[[11]], l[[12]], l[[13]], l[[14]], l[[15]],
          l[[16]], l[[17]], l[[18]], l[[19]], l[[20]],
          l[[21]], nrow=5, ncol=5
          )

ggsave('~/Desktop/plot.pdf', plot=p,  width = 30, height = 30, units = c("in"))


