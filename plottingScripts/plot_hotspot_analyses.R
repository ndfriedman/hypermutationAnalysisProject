#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

plot_hotspot_utilization_curve <- function(df, title, hideLegend = TRUE){
  plt <- ggplot(df, aes(x=percentile, y=val, group=class))+
    geom_path(aes(colour=class))+
    #geom_segment(x=0, xend=1, y=0, yend=1, alpha=0.1, colour='gray')+
    ylab('Fraction of all hotspot mutations')+
    xlab('Common                 Rare')+
    xlim(0,1)+
    ylim(0,1)+
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank())+
    labs(colour = "Gene Type")+
    ggtitle(title)
  if(hideLegend){
    plt <- plt + theme(legend.position = "none")
  }
  return(plt)
}

#TEMP TEMP TEMP
dfEndo <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/EndoHotspotUtilizationCurve.tsv',sep = '\t', header=TRUE)
dfColo <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/ColoHotspotUtilizationCurve.tsv',sep = '\t', header=TRUE)
dfGlio <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/GlioHotspotUtilizationCurve.tsv',sep = '\t', header=TRUE)

arrangedPlot <- plot_grid(
  plot_hotspot_utilization_curve(dfEndo, 'Endometrial Cancer'),
  plot_hotspot_utilization_curve(dfColo, 'Colorectal Cancer'),
  plot_hotspot_utilization_curve(dfGlio, 'Glioma'),
  cowplot::get_legend(plot_hotspot_utilization_curve(dfGlio, 'Glioma', hideLegend = FALSE)), #THIS IS THE LEGEND ONLY
  ncol = 4, rel_widths = c(1,1,1,.5)
)
fullPlot <- plot_grid(
                      ggplot() + ggtitle('Hotspot Utilization in Hypermutated Cancer'),
                      arrangedPlot,
                      ggplot() + labs(caption = 'plot_hotspot_analyses.R  hotspot_utilization_analysis.ipynb'),
                      nrow=3, rel_heights = c(.1,1,.1))

ggsave('~/Desktop/plot.pdf', plot=fullPlot,  width = 10, height = 4, units = c("in"))


plot_log_odds_p_info <- function(df, title, plotColour){
  
  qobj <- qvalue(p = df$pvalueSIGMuts, pi0 = 1)
  df$qvalue <- qobj$qvalues
  
  plt <- ggplot(df)+
    #geom_point(aes(x=-1*log(qvalue), y=log(oddsratioSIGMuts)), colour=plotColour, alpha=.75)+
    geom_point(aes(x=-1*log(pvalueSIGMuts), y=log(oddsratioSIGMuts)), colour=plotColour, alpha=.75)+
    
    #geom_point(aes(x=-1*log(pvalueSIGMuts), y=log(oddsratioSIGMuts), colour=mutOfInterestFrac), alpha=.75)+
    #scale_colour_viridis_c()+
    
    geom_vline(aes(xintercept= -1*log(.05)))+
    geom_hline(aes(yintercept= log(1)))+
    ylab('log oddsratio')+
    xlab('-log q value')+
    ylim(0,6)+
    xlim(0,20)+
    ggtitle(title)
  return(plt)
}

dfHotspotPOLE <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/hotspotPValues_pole.tsv',sep = '\t', header=TRUE)
hotspotPlotPOLE <- plot_log_odds_p_info(dfHotspotPOLE, 'Hotspot muts', 'green')
dfOncPOLE <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/oncogenicPValues_pole.tsv',sep = '\t', header=TRUE)
oncPlotPOLE <- plot_log_odds_p_info(dfOncPOLE, 'Oncogenic muts', 'green')

dfHotspotTMZ <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/hotspotPValues_tmz.tsv',sep = '\t', header=TRUE)
hotspotPlotTMZ <- plot_log_odds_p_info(dfHotspotTMZ, 'Hotspot muts', 'blue')
dfOncTMZ <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/oncogenicPValues_tmz.tsv',sep = '\t', header=TRUE)
oncPlotTMZ <- plot_log_odds_p_info(dfOncTMZ, 'Oncogenic muts', 'blue')

dfHotspotAPOBEC <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/hotspotPValues_apobec.tsv',sep = '\t', header=TRUE)
hotspotPlotAPOBEC <- plot_log_odds_p_info(dfHotspotAPOBEC, 'Hotspot muts', 'red')
dfOncAPOBEC <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/oncogenicPValues_apobec.tsv',sep = '\t', header=TRUE)
oncPlotAPOBEC <- plot_log_odds_p_info(dfOncAPOBEC, 'Oncogenic muts', 'red')

alignedPlotPOLE <- plot_grid(hotspotPlotPOLE, oncPlotPOLE)
alignedPlotTMZ <- plot_grid(hotspotPlotTMZ, oncPlotTMZ)
alignedPlotAPOBEC <- plot_grid(hotspotPlotAPOBEC, oncPlotAPOBEC)

fullPlot <- 
    plot_grid(ggplot() + ggtitle('Sig specific motifs and enrichment for activating mutations'),
    ggplot()+ggtitle('POLE Hypermutated Cases (n=77)'),
    alignedPlotPOLE,
    ggplot()+ggtitle('TMZ Hypermutated Cases (n=59)'),
    alignedPlotTMZ,
    ggplot()+ggtitle('APOBEC High mut burden: >75mut_mb (n=93)'),
    alignedPlotAPOBEC,
    ggplot()+labs(caption = 'plot_hotspot_analyses.R, hotspot_and_motif_analysis.ipynb'),
    nrow=8, rel_heights = c(.1,.1, 1, .1, 1,.1, 1, .1)
)


ggsave('~/Desktop/plot.pdf', plot=fullPlot,  width = 8, height = 10, units = c("in"))


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pvalues <- dfOncPOLE$pvalueSIGMuts
qobj <- qvalue(p = pvalues, pi0 = 1)
plot(qobj)


print(pvalues - qobj$qvalues)
#####################################

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/clonalVsSubclonalOncMutUtilization_oncogenic.tsv',sep = '\t', header=TRUE)
p1 <- ggplot(df, aes(x=rank, y=frac, group=class, colour=class))+
  geom_line()+
  ggtitle('Oncogenic muts')

df2 <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/clonalVsSubclonalOncMutUtilization_hotspot.tsv',sep = '\t', header=TRUE)
p2 <- ggplot(df2, aes(x=rank, y=frac, group=class, colour=class))+
  geom_line()+
  ggtitle('Hotspot muts')

alignedP <- plot_grid(p1, p2, ncol=2)
fullPlot <- plot_grid(ggplot()+ggtitle('Gene utilization Clonal vs Subclonal'),
  alignedP, ggplot() + labs(caption='gene_multiple_hit_analysis.ipynb, plot_hotspot_analyses.R'), nrow=3, rel_heights = c(.1,1,.1))

ggsave('~/Desktop/plot.pdf', plot=fullPlot,  width = 10, height = 4, units = c("in"))



###########
##############
#################
#####################
##########################

make_fractional_panels <-function(df, brks, legendTitle, panelLabel){
  df$group <- cut(df$Nmut, breaks=brks)
  
  newDf <- data.frame(matrix(ncol = 4, nrow = 0))
  colNames <- c("TMB_group", "SigMotifRelatedMean", "NotSigMotifRelatedMean", "Nmut_Mean")
  names(newDf) <- colNames
  for(g in unique(df$group)){
    groupMeanRelated <- mean(df[df$group == g,]$fracRelatedMotif, na.rm=TRUE)
    groupMeanUnrelated <- mean(df[df$group == g,]$fracRelatedNotMotif, na.rm=TRUE)
    nmutMean <- mean(df[df$group == g,]$Nmut, na.rm=TRUE)
    localDf <- data.frame(g, groupMeanRelated, groupMeanUnrelated, nmutMean)
    names(localDf) <- colNames
    newDf <- rbind(newDf, localDf)
  }
  
  panelPlot <- ggplot(newDf)+
    geom_tile(aes(x=reorder(TMB_group, Nmut_Mean), y = .5, fill=SigMotifRelatedMean))+
    geom_tile(aes(x=reorder(TMB_group, Nmut_Mean), y = -.5, fill=NotSigMotifRelatedMean))+
    scale_fill_viridis_c(option='magma', direction=-1)+
    theme(axis.line = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank(),
          axis.text.y = element_blank(),axis.ticks.y = element_blank(), axis.title.y=element_text(size=6),
    )+
    ylab(panelLabel)+
    guides(fill=guide_legend(title=legendTitle))
  
  return(panelPlot)
}

plot_motif_analysis_with_panels <- function(df, brksm, title='Hotspot mutations in non-MSI Endometrial Cancer', 
                                            motifName='Pole motif\n mutations', 
                                            notMotifName='Not pole motif\n mutations',
                                            legendTitle='Frac Hotspots\n in Endo-Cancer\n genes',
                                            panelLabel='Other motifs  POLE motif'){
  chart <- ggplot(df)+
    
    stat_summary_bin(aes(x=Nmut, y=hotspotMutsNotAtSigMotif, colour=notMotifName),
                     breaks=brks, alpha=0.75)+
    stat_summary_bin(aes(x=Nmut, y=nHotspotsAtMotif, colour=motifName),
                     breaks=brks, alpha=0.75)+
    scale_x_log10(breaks=brks)+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+
    ylab('N hotspot mutations at Sig Motif')+
    ylim(0,40)+
    guides(fill=guide_legend(title='Mutation context'))
  
  panel <- make_fractional_panels(df, brks, legendTitle, panelLabel)
  alignedPlot <- plot_grid(
    ggplot()+ggtitle(title),
    chart, panel,
    ggplot()+labs(caption='hotspot_and_motif_analysis.ipynb   plot_hotspot_analyses.R'),
    nrow=4, rel_heights = c(.1,1,.2,.1))
  return(alignedPlot)
}

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/poleEndometrialHotspotCharacteristics.tsv',sep = '\t', header=TRUE)
brks = c(0,5,10,20,40,75,125,200,300,max(df$Nmut))
plt <- plot_motif_analysis_with_panels(df,brks)

dfBladder <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/apobecBladderHotspotCharacteristics.tsv',sep = '\t', header=TRUE)
dfGlioma <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/tmzGliomaHotspotCharacteristics.tsv',sep = '\t', header=TRUE)

bladderBrks = c(0,5,10,15,30,45,75,100,max(dfBladder$Nmut))
plot_motif_analysis_with_panels(dfBladder, bladderBrks,
  title='Hotspot mutations in Bladder Cancer', 
  motifName='APOBEC motif\n mutations', 
  notMotifName='Not APOBEC motif\n mutations',
  legendTitle='Frac Hotspots\n in Bladder-Cancer\n genes')


gliomaBrks = c(0,5,10,20,40,75,125,200,max(dfGlioma$Nmut))
p <- plot_motif_analysis_with_panels(dfGlioma, gliomaBrks,
                                title='Hotspot mutations in Glioma', 
                                motifName='TMZ motif\n mutations', 
                                notMotifName='Not TMZ motif\n mutations',
                                legendTitle='Frac Hotspots\n in Glioma\n genes',
                                panelLabel='Other motifs  TMZ motif')


#####
############
###############
####################
#########################
##############################
######################################
#############################################
######################################
################################
##########################
###################
###############
###########
####
emptyTheme <- theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())

plot_mutation_commonality <- function(df, cancerType, textVal, yc = .8){
    #ggplot(df[df$cancer_type == cancerType,], aes(x=mutPercentile, y=1-ranking))+
    #stat_summary_bin(bins=15)+
    ggplot(df[df$cancer_type == cancerType,], aes(x=mutPercentile, y=1-ranking, colour=mutNumber))+
      geom_point()+
  
    #geom_text(aes(x=0.5, y=yc, label=textVal), size=2.5)+
    ylim(0,1)+
    xlab('<-----most common onc mutation in case\nleast common onc mutation in case----->')+
    ylab('<---------uncommon genes               common genes------------>')+
    theme(legend.position = 'None')+
    theme(axis.ticks = element_blank(), axis.text = element_blank())+
    geom_segment(aes(x=0, y=1, xend=1, yend=0), colour='black')+
    ggtitle(cancerType)+
    emptyTheme
}

plot_mutation_commonality_vus <- function(df, cancerType){
  ggplot(df[df$cancer_type == cancerType,], aes(x=mutPercentile, y=1-ranking))+
    stat_summary_bin(bins=15)+
    ylim(0,1)+
    xlab('<----- ----->')+
    ylab('<--------- ------------>')+
    theme(axis.ticks = element_blank(), axis.text = element_blank())+
    geom_segment(aes(x=0, y=1, xend=1, yend=0), colour='black')+
    ggtitle(cancerType)+
    emptyTheme
}

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/plotOncogenicOccurrenceFreq.tsv',sep = '\t', header=TRUE)
dfVUS <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/plotVUSOccurrenceFreq.tsv',sep = '\t', header=TRUE)
dfHotspot <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/plotHOTSPOTOccurrenceFreq.tsv',sep = '\t', header=TRUE)


endoLabel <- paste(median(df[df$cancer_type == 'Endometrial Cancer',]$nConditionMuts)/2, 'th mutation in median case')
endoPlot <- plot_mutation_commonality(df, 'Endometrial Cancer', endoLabel)
coloLabel <- paste(median(df[df$cancer_type == 'Colorectal Cancer',]$nConditionMuts)/2, 'th mutation in median case')
coloPlot <- plot_mutation_commonality(df, 'Colorectal Cancer', coloLabel, yc=.775)
gliomaLabel <- paste(median(df[df$cancer_type == 'Glioma',]$nConditionMuts)/2, 'th mutation in median case')
gliomaPlot <- plot_mutation_commonality(df, 'Glioma', gliomaLabel)

alignedPlot <- plot_grid(endoPlot, coloPlot, gliomaPlot, ncol=3)
leg <- get_legend(ggplot(df, aes(x=mutPercentile, y=1-ranking, group=nmutCategory, colour=nmutCategory))+
  stat_summary())
  
#alignedPlotWithLegend <- plot_grid(alignedPlot, leg, ncol=2, rel_widths = c(1,.2))
alignedPlotWithLegendAndTitle <- plot_grid(ggplot() + ggtitle('Oncogenic mutations in hypermutated cancer') +theme(plot.title = element_text(hjust=.5, face='bold')),
                                           alignedPlot, ggplot() + labs(caption='plot_hotspot_analyses.R  hotspot_utilization_analysis.ipynb'), nrow=3, rel_heights = c(.05,1, .025))



ggsave('~/Desktop/plot.pdf', plot=alignedPlotWithLegendAndTitle,  width = 12, height = 5.5, units = c("in"))


#EXPLANATORY PLOT
oncPlot <- ggplot(df[(df$cancer_type == 'Endometrial Cancer') & (df$caseNmut > 350),], aes(x=mutNumber, y=1-ranking, group=Tumor_Sample_Barcode, colour=Tumor_Sample_Barcode))+
  geom_path(alpha=0.5)+
  ylim(0,1)+
  ylab('<---- uncommon genes | common genes ------>')+
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())+
  emptyTheme+
  theme(legend.position = "none")+
  ggtitle('Oncogenic mutations')

vusPlot <- ggplot(dfVUS[(dfVUS$cancer_type == 'Endometrial Cancer') & (dfVUS$caseNmut > 350),], aes(x=mutNumber, y=1-ranking, group=Tumor_Sample_Barcode, colour=Tumor_Sample_Barcode))+
  geom_path(alpha=0.5)+
  ylim(0,1)+
  ylab('<---- uncommon genes | common genes ------>')+
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())+
  emptyTheme+
  theme(legend.position = "none")+
  ggtitle('VUS mutations')

hotspotPlot <- ggplot(dfHotspot[(dfHotspot$cancer_type == 'Endometrial Cancer') & (dfHotspot$caseNmut > 350),], aes(x=mutNumber, y=1-ranking, group=Tumor_Sample_Barcode, colour=Tumor_Sample_Barcode))+
  geom_path(alpha=0.5)+
  ylim(0,1)+
  ylab('<---- uncommon genes | common genes ------>')+
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())+
  emptyTheme+
  theme(legend.position = "none")+
  ggtitle('Hotspot mutations')

alignedPlot <- plot_grid(vusPlot, oncPlot, hotspotPlot, ncol=3)
fullPlot <- plot_grid(ggplot() + ggtitle('mutation recurrence by mutation number (MOST ULTRA MUTATED ENDOMETRIAL)'), 
                                         alignedPlot, ggplot() + labs(caption='plot_hotspot_analyses.R'), nrow=3, rel_heights=c(.1,1,.05))

ggsave('~/Desktop/plot.pdf', plot=fullPlot,  width = 8, height = 4.5, units = c("in"))

