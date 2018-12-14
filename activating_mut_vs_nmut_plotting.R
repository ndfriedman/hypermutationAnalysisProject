#written by noah friedman
#a template for R scripts for plotting
#A script that plots the pan cohort numbers of activating mutations by mutation burden and the ratio of activating mutations to total number of mutations

library(ggplot2)
library(grid)
require(cowplot)
library(egg)


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


if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

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


