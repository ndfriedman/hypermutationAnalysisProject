#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)
library(ggrepel)


if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

break_axis <- function(y, maxlower, minupper=NA, lowerticksize, upperticksize, ratio_lower_to_upper) {
  if(is.na(minupper)) {
    breakpos <- maxlower
    lowerticklabels <- seq(0,breakpos,by=lowerticksize); lowerticklabels
    upperticklabels <- seq(breakpos+upperticksize,max(y)+upperticksize,by=upperticksize); upperticklabels
    ticklabels <- c(lowerticklabels, upperticklabels); ticklabels
    lowertickpos <- lowerticklabels
    uppertickspacing <- ratio_lower_to_upper * lowerticksize
    uppertickpos <- breakpos + ((1:length(upperticklabels))*uppertickspacing)
    tickpos <- c(lowertickpos, uppertickpos)
    newy <- as.numeric(y)
    ind <- newy > breakpos
    newy[ind] <- breakpos + uppertickspacing*((newy[ind]-breakpos) / upperticksize)
    list(newy=newy, breaks=tickpos, labels=ticklabels, limits=range(tickpos))
  } else {
    lowerticklabels <- seq(0,maxlower,by=lowerticksize); lowerticklabels
    upperticklabels <- seq(minupper,max(y)+upperticksize,by=upperticksize); upperticklabels
    ticklabels <- c(lowerticklabels, upperticklabels); ticklabels
    lowertickpos <- lowerticklabels
    uppertickspacing <- ratio_lower_to_upper * lowerticksize
    uppertickpos <- maxlower + 0.5*lowerticksize + ((1:length(upperticklabels))*uppertickspacing)
    tickpos <- c(lowertickpos, uppertickpos)
    newy <- as.numeric(y)
    ind <- newy > maxlower
    newy[ind] <- maxlower + 0.5*lowerticksize + 1*uppertickspacing + uppertickspacing*((newy[ind]-minupper) / upperticksize)
    list(newy=newy, breaks=tickpos, labels=ticklabels, limits=range(tickpos))
  }
}


emptyTheme <- theme(axis.line = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())

df <- read.table('/Users/friedman/Desktop/WORK/msiSiteOccurenceInfo.tsv',sep = '\t', header=TRUE)

df[(df$repeat_length == 4) & (df$nOccurences > 15),]

#MAKE THE BREAKS
x <- break_axis(df$nOccurences, maxlower=175, minupper=300, lowerticksize=25, upperticksize=25, ratio_lower_to_upper=0.5)
df$newBreaks <- x$newy

plt <- ggplot(df[(df$isTumorSuppresor == "True") | (df$isTumorSuppresor == "False"),], aes(x=repeatLengthsLabel, y=newBreaks))+
  geom_boxplot(fatten = NULL, outlier.shape = NA)+ #TODO FIX THIS TO PLOT THE MIDDLE LINE AS THE MEAN!!!
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 0.75, size = 1, linetype = "solid")+
  #geom_jitter(aes(colour=isTumorSuppresor), shape=16, position=position_jitter(0.1), alpha=0.5)+
  geom_jitter(aes(colour=basePair), shape=16, position=position_jitter(0.1), alpha=0.75)+
  
  #ALEX axis breaks
  scale_y_continuous(breaks=x$breaks, labels=x$labels, limits=x$limits)+
  geom_text_repel(aes(label=label),
                  size=2,
                  segment.size  = 0.5,
                  segment.color = "grey50",
                  direction     = "both")+
  ggtitle('Frequency of MSI Indels by Repeat Length')+
  emptyTheme+
  xlab('Repeat Length')+
  ylab('N cases with mutation in IMPACT')+
  theme(plot.title = element_text(face = 'bold'))

plotWithCaption <- plot_grid(plt, ggplot()+labs(caption='msi_mut_features_analysis.ipynb, plotMSIAnalyses.R'),
                             nrow=2, rel_heights = c(1,.05))
ggsave('~/Desktop/plot.pdf', plot=plotWithCaption,  width = 5, height = 5, units = c("in"))
 
#
#
#

#
########
############
########3######
plt <- ggplot(df[(df$repeat_length > 4) & (df$repeat_length < 9) & (df$basePair != '_other'),], aes(x=reorder(basePairAndLength, orderingVal), y=nOccurences))+
  stat_summary(aes(colour=basePair))+
  ggtitle('Frequency of MSI Indels by Repeat Length')+
  emptyTheme+
  xlab('Repeat Length and Mutated Base')+
  ylab('N cases with mutation in IMPACT')+
  theme(axis.text.x = element_text(angle=70, hjust=1))
  #coord_cartesian(ylim=c(0,50))

ggsave('~/Desktop/plot.pdf', plot=plt,  width = 5, height = 5, units = c("in"))

#
####
#########
###############
#################
#######################
##############################

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/doubleMutationMSIInfo.tsv', sep = '\t', header=TRUE)

#alex g break axis
x <- break_axis(df$nOccurences, maxlower=15, minupper=42, lowerticksize=1, upperticksize=2, ratio_lower_to_upper=0.5)
df$x <- x$newy

p <- ggplot(df, aes(x=repeat_length, y=df$x))+
  geom_jitter(aes(colour=basePair), size=1, position=position_jitter(0.1))+
  emptyTheme+
  geom_text_repel(aes(label=label),
                  size=4,
                  segment.size  = 0.5,
                  segment.color = "grey50",
                  direction     = "both")+
  scale_y_continuous(breaks=x$breaks, labels=x$labels, limits=x$limits)+
  geom_hline(aes(yintercept=15), linetype=12345678)+
  ylab('N cases')+
  ggtitle('Double Allele Hit Indel Muts')+
  theme(plot.title = element_text(face = 'bold'))

plotWithCaption <- plot_grid(p, ggplot()+labs(caption='msi_mut_features_analysis.ipynb, plotMSIAnalyses.R'),
                             nrow=2, rel_heights = c(1,.05))
ggsave('~/Desktop/plot.pdf', plot=plotWithCaption,  width = 5, height = 5, units = c("in"))


#
########
#############
###################
#########################
##############################
####################################
#########################################

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/msiIndelFeatures.tsv', sep = '\t', header=TRUE)
basePairPlot <- ggplot(df)+
  geom_histogram(aes(x=cg_atRATIO, fill=dominantSig), bins=15)+
  scale_fill_brewer(palette = "Set2")+
  emptyTheme+
  scale_x_log10(breaks=c(.5,1,1.5,2,3,5))+
  ylab('n cases')+
  theme(legend.position = "none")+
  xlab('Ratio of C-G to A-T MSI sites targeted')+
  ggtitle('Base Pair Composition of Cases')+
  geom_vline(aes(xintercept=1), linetype=12345678)+
  theme(axis.text.x = element_text(angle=70, hjust=1))

basePairGenePlot <- ggplot(df[df$mmrGene != 'None',])+
  geom_histogram(aes(x=cg_atRATIO, fill=mmrGene), bins=10)+
  scale_fill_brewer(palette = "Set2")+
  emptyTheme+
  scale_x_log10(breaks=c(.5,1,1.5,2,3,5))+
  ylab('n cases')+
  #theme(legend.position = "none")+
  xlab('Ratio of C-G to A-T MSI sites targeted')+
  ggtitle('Base Pair Composition of Cases by Gene')+
  geom_vline(aes(xintercept=1), linetype=12345678)+
  theme(axis.text.x = element_text(angle=70, hjust=1))

indelLengthPlot <-  ggplot(df)+
    geom_histogram(aes(x=averageIndelLen, fill=dominantSig), bins=20)+
    scale_fill_brewer(palette = "Set2")+
    emptyTheme+
    ylab('n cases')+
    guides(fill=guide_legend(title='Dominant Signature'))+
    xlab('Average repeat length of MSI Indel')+
    ggtitle('Repeat Length Characteristics of Cases')+
  theme(axis.text.x = element_text(angle=70, hjust=1))

alignedPlot <- plot_grid(basePairPlot, indelLengthPlot, ncol=2)
fullPlot <- plot_grid(alignedPlot, ggplot()+labs(caption='msi_mut_features_analysis.ipynb plotMSIAnalyses.R'),
                      nrow=2, rel_heights = c(1,.05))

ggsave('~/Desktop/plot.pdf', plot=basePairGenePlot,  width = 3.5, height = 3.5, units = c("in"))



#
###
#######
###########
###############
##################
######################
##################
###############
#############
##########
########
####
#

emptyTheme <- theme(axis.line = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/frameShiftTestMsi.tsv', sep = '\t', header=TRUE)

p <- ggplot(df, aes(x=-1*log(pValue), y=log(oddsRatio)))+
  geom_point(aes(colour=isTumorS))+
  geom_text_repel(aes(label = displayLabel))+
  #scale_x_log10()+
  #scale_y_log10()+
  emptyTheme+
  ggtitle('Fisher test: MSI frameshift vs not frameshift obs/expected')+
  geom_hline(aes(yintercept=log(1)))+
  ylab('Log odds ratio')+
  xlab('Log p value')

ggsave('~/Desktop/plot.pdf', plot=p,  width = 6, height = 5, units = c("in"))


#
##
######
############
######################
#############################
#####################################
############################################
#######################################
#############################
########################
###############
########
###
#
ggplot(df, aes(x=reorder(basePair, orderingVal), y=nOccurences))+
  stat_summary(aes(colour=basePair))



df <- read.table('/Users/friedman/Desktop/WORK/msiSiteOccurenceInfo.tsv', sep = '\t', header=TRUE)

df$isTumorSuppresor

ggplot(df[(df$basePairClass != 'other') & (df$repeat_length < 9),], aes(x=reorder(class, orderingVal),
        y=nOccurences, group=isTumorSuppresor, colour=isTumorSuppresor))+
  stat_summary(geom="line", fun.y="mean")+
  stat_summary()+
  scale_y_continuous(breaks=c(0,10,20,30,40,50,75))+
  theme(axis.text=element_text(angle=70, hjust=1))+
  ylab('Average Occurences per Allele')+
  xlab('Repeat length and basepair')+
  ggtitle('MSI Indels preferentially mutate\ntumor suppressors')+
  labs(caption='plotMSIAnalyses.R\ngenerate_we_see_more_mutations_than_expected_figure.ipynb')
  

df <- read.table('/Users/friedman/Desktop/WORK/poleSitOccurenceInfo.tsv', sep = '\t', header=TRUE)
ggplot(df, aes(x=reorder(penta, orderingVal),y=nmut, group=tumorSup, color=tumorSup))+
  stat_summary()+
  theme(axis.text=element_text(angle=70, hjust=1))+
  xlab('Pentanucleotide Context')+
  ylab('Average truncating mutations per allele\n across cohort')+
  labs(caption='plotMSIAnalyses.R\nallele_occurences_at_signature_favored_motifs.ipynb')+
  ggtitle('POLE signature induced truncating mutations\npreferentially target tumor suppressors')

  
  
  
  
  
  
  
  
  
  










