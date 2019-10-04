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

emptyTheme <- theme(axis.line = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())


df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/alleleSpecificity.tsv',sep = '\t', header=TRUE)
ggplot(df[(df$nSteps > 0) | (df$diffScore > 0),], aes(nSteps, y=diffScore))+
  geom_point(aes(size=nMutations), colour='orange')+
  geom_text_repel(aes(label=displayName))+
  geom_text(aes(x = max(df$nSteps), y= max(df$diffScore)), 
            label='only one allele:\nTP53(n=18)\nSMAD4(n=7)\nPARK2(n=2)\nLATS2(n=2)')+
  xlim(0,20)+
  ylim(0,1.75)+
  xlab('Allele Specificity')+
  ylab('Distribution Divegence From Expectation')+
  ggtitle('Allele Specificity of \nTruncating TS-Mutations at\n POLE Favored Sites')+
  emptyTheme+
  labs(caption='plot_allele_mutation_information.R\n allele_occurences_at_signature_favored_motifs.ipynb')

ggsave('~/Desktop/plot.pdf', plot=plt,  width = 1, height = 1, units = c("in"))



#
#
#######
############
###########################
########################
############
###
#

plot_penta_alleles <- function(df, penta){
  df <- df[df$penta == penta,]
  yVal <- dim(df[df$count == 0,])[1]/3
  plt <- ggplot(df, aes(x=count))+
    geom_bar(stat='count')+
    geom_text_repel(aes(label=displayLabel, y=0), nudge_y=yVal, force=10)+
    xlab('N occurences of allele in cohort')+
    ggtitle(penta)
  return(plt)
}

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/alleleCounts.tsv',sep = '\t', header=TRUE)

p1 <- plot_penta_alleles(df, 'TT(C>T)GA')
p2 <- plot_penta_alleles(df, 'TT(C>A)TC')
p3 <- plot_penta_alleles(df, 'TT(C>A)TG')
p4 <- plot_penta_alleles(df, 'AT(C>T)GA')
p5 <- plot_penta_alleles(df, 'GT(C>T)GA')
p6 <- plot_penta_alleles(df, 'TT(C>A)TT')

fullPlotPentaTrunc <- plot_grid(p1, p2, p3, p4, p5, p6, nrow=3, ncol=3)
fullPlotWithTitle <- plot_grid(ggplot()+ggtitle('Recurrence of truncating mutations across favored POLE penta-nucleotides in all POLE cases'),
                               fullPlotPentaTrunc,
                               ggplot()+ labs(caption='plot_allele_mutation_information.R\nsummarize_mutations_seen_and_never_seen.ipynb'),
                              nrow=3, rel_heights = c(.1, 1, .05))

ggsave('~/Desktop/plot.pdf', plot=fullPlotWithTitle,  width = 10, height = 10, units = c("in"))


#
df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/hotspotAlleleCounts.tsv',sep = '\t', header=TRUE)
p1 <- plot_penta_alleles(df, 'TT(C>A)TT')
p2 <- plot_penta_alleles(df, 'TT(C>T)GA')
p3 <- plot_penta_alleles(df, 'TT(C>A)TC')
p4 <- plot_penta_alleles(df, 'TT(C>T)GG')
p5 <- plot_penta_alleles(df, 'TT(C>T)GT')
p6 <- plot_penta_alleles(df, 'TT(C>T)GC')
p7 <- plot_penta_alleles(df, 'TT(C>A)TG')
p8 <- plot_penta_alleles(df, 'AT(C>A)TT')
p9 <- plot_penta_alleles(df, 'AT(C>T)GA')
p10 <- plot_penta_alleles(df, 'GT(C>T)GA')

fullPlotPentaTrunc <- plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, nrow=4, ncol=3)
fullPlotWithTitle <- plot_grid(ggplot()+ggtitle('Recurrence of hotspot mutations across favored POLE penta-nucleotides in all POLE cases'),
                               fullPlotPentaTrunc, ggplot()+ labs(caption='plot_allele_mutation_information.R\nsummarize_mutations_seen_and_never_seen.ipynb'),
                               nrow=3, rel_heights = c(.1, 1. ,.05))

ggsave('~/Desktop/plot.pdf', plot=fullPlotWithTitle,  width = 10, height = 10, units = c("in"))



#
##
####
#######
#######
############
######
######
#
#

#PLOT double mutation as it relates to stop gain mutation alleles
df2 <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/stopGainAlleleRecurrenceAndDoubles.tsv',sep = '\t', header=TRUE)
ggplot(df2, aes(x= reorder(displayAllele, orderingVal), y=multipletPresentInCase))+
  stat_summary(aes(colour=geneClass))+
  #geom_point(aes(colour = geneClass))+
  theme(axis.text.x = element_text(angle=90))+
  ylab('Fraction of cases harboring allele\nwith additional mutation in gene')+
  xlab('allele')+
  ggtitle('co-occurence of POLE favored stop-gain alleles\n and second hit mutations in the same gene')+
  labs(caption='summarize_mutations_seen_and_never_seen.ipynb\nplot_allele_mutation_information.R')

df2$displayAllele
#ggplot(df, aes(x=fullClass, y=frac))+
#  geom_boxplot()

df$orderingVal

