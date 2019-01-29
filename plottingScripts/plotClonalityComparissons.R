#written by noah friedman
#a template for R scripts for plotting
#Does simple comparisson plots of hte clonality fraction by mutation types

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(stringr)
library(ggpubr)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")


#a little function that allow us to label bars with n occurences
n_fun <- function(x){
  adjust <- 0.05
  return(data.frame(y = mean(x) + adjust, label = paste0("n = ",length(x))))
}

do_plot_details <- function(df, xValParam, xSortParam, colors, title=""){
  df <- my.rename(df, xValParam, "xVal")
  df <- my.rename(df, xSortParam, "xSort")
  
  ggplot(df)+
  stat_summary_bin(data=df, aes(x=xVal, y=clonal, fill=xSort), fun.y = "mean", geom = "bar",  position="dodge")+
  #stat_summary(aes(x = xVal, y = clonal, fill=xSort), fun.data = n_fun, geom = "text")+
  scale_fill_manual(values=colors)+
  ylab("fraction of clonal muts")+
  ylim(0,1)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  xlab("gene mutated")+
  ggtitle("Clonal fraction of muts in TMZ by motif and gene")
}

do_simple_plot <- function(df, xValParam, title="", xlabel="", compYPos=NA){
  df <- my.rename(df, xValParam, "xVal")
  plt <- ggplot(df)+
    stat_summary(data=df, aes(x = xVal, y = clonal), fun.y = "mean", geom = "bar")+
    stat_summary(aes(x = xVal, y = clonal), fun.data = n_fun, geom = "text")+
    ggtitle(title)+
    xlab(xlabel)+
    ylab("fraction of clonal muts")+
    ylim(0,1)+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 15))
  
  stat.test <- compare_means(
    clonal ~ xVal, data = df,
    method = "t.test"
  )
  stat.test <- stat.test %>%
    mutate(y.position = compYPos) #we pass in the positions for the comparissons here
  plt <- plt + stat_pvalue_manual(stat.test, label = "p.signif")
  
  return(plt)
}

#########GLIOMA PLOTTING
gliomaDf <- read.table('~/Desktop/dataForLocalPlotting/gliomaClonalityData.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
tmzColors = c("pink", "#ADD8E6", "#800000", "#FF0000", "Black")
#ORDER them and use the exact names we have
gliomaDf$geneClass = factor(gliomaDf$geneClass, levels = c("TERT-TP53-IDH1-PTEN-EGFR", "ATRX-NF1-CIC-PIK3CA-PIK3CR1-RB1-NOTCH1-PDGFRA-BRAF-PTPN11", "other genes")) #properly order the charts

#DO a couple different singleton plots

pltAllMutsByClass <- do_simple_plot(gliomaDf, "geneClass", title = "Fraction of clonal muts in genes TMZ-glioma")
pltHotspotMutsByClass <- do_simple_plot(gliomaDf[gliomaDf$mutEffectLabel == 'hotspot',], "geneClass", title = "Fraction of clonal hotspot muts TMZ-glioma")

pos <- c(.7, .75,.8,.85,.9,.95,1,1.05,1.1, 1.15)
pltAllMutsByTMZMotif <- do_simple_plot(gliomaDf[gliomaDf$mutEffectLabel == 'hotspot',], "sigMotifLabel", title = "Fraction of clonal muts by TMZ motif", compYPos=pos)
pos <- c(.8,.9,1)
pltTMZTopMutsByGenes <- do_simple_plot(gliomaDf[gliomaDf$sigMotifLabel == 'topTMZ',], "geneClass", title = "Clonality in genes mutated in TMZ most favored motif ", compYPos=pos)
pltAgingTopMutsByGenes <- do_simple_plot(gliomaDf[gliomaDf$sigMotifLabel == 'fourthTMZ',], "geneClass", title = "Clonality in genes mutated in Aging most favored motif ", compYPos=pos)

fullPlot <- do_plot_details(gliomaDf, "geneClass", "sigMotifLabel", tmzColors, title="")

################################COLORECTAL

colorectalDf <- read.table('~/Desktop/dataForLocalPlotting/colorectalClonalityData.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
pltAllMutsByClassColorectal <- do_simple_plot(colorectalDf, "geneClass", title = "Fraction of clonal muts Hypermuated Colorectal", compYPos=pos)
pos <- c(.8)
pltAllMutsByMotifColorectal <- do_simple_plot(colorectalDf, "sigMotifLabel", title = "Fraction of clonal muts Hypermuated Colorectal", xlabel="Motif of trinucleotide mutation", compYPos=pos)

#######################################
#################MELANOMA
melanomaDf <- read.table('~/Desktop/dataForLocalPlotting/melanomaClonalityData.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
pos <- c(.8,.9,1)
pltAllMutsByClassMelanoma <- do_simple_plot(melanomaDf, "geneClass", title = "Fraction of clonal muts Hypermuated Melanoma",xlabel="", compYPos=pos)
pos <- c(.8)
pltAllMutsByClassMelanoma <- do_simple_plot(melanomaDf, "sigMotifLabel", title = "Fraction of clonal muts Hypermuated Melanoma",xlabel="Motif of trinucleotide mutation", compYPos=pos)

###########################################
##########################LUNG CANCER
lungDf <- read.table('~/Desktop/dataForLocalPlotting/lungClonalityData.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
pos <- c(.8,.9,1)
pltAllMutsByClassLung <- do_simple_plot(lungDf, "geneClass", title = "Fraction of clonal muts high mut burden Lung Cancer",xlabel="Motif of trinucleotide mutation", compYPos=pos)
pos<-c(.8)
pltAllMutsByClassLung <- do_simple_plot(lungDf, "sigMotifLabel", title = "Fraction of clonal muts high mut burden Lung Cancer",xlabel="Motif of trinucleotide mutation", compYPos=pos)

######################################
################BLADER CANCER
bladderDf <- read.table('~/Desktop/dataForLocalPlotting/bladderClonalityData.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
pos <- c(.8,.9,1)
pltAllMutsByClassBladder <- do_simple_plot(bladderDf, "geneClass", title = "Fraction of clonal muts high mut burden Bladder Cancer",xlabel="Motif of trinucleotide mutation", compYPos=pos)
pos<-c(.8)
pltAllMutsByClassBladder <- do_simple_plot(bladderDf, "sigMotifLabel", title = "Fraction of clonal muts high mut burden Bladder Cancer",xlabel="Motif of trinucleotide mutation", compYPos=pos)







ggsave('~/Desktop/testNoah1.pdf', plot=plt1,  width = 8, height = 5, units = c("in"))

