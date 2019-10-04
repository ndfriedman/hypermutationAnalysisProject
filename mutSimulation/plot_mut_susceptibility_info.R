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


df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/hotspotObsExpectedPOLE.tsv',sep = '\t', header=TRUE)

ggplot()+
  geom_path(data=df, aes(x=nmutIM341, y=nHotspotsExpected, colour='Hotspots Expected\nmethod 2'))+
  geom_smooth(data=dfPenta, aes(x=nmut, y=nExpectedHotspots, colour='Hotspots Expected\nmethod 1'))+
  geom_smooth(data=df, aes(x=nmutIM341, y=nHotspotsInducedByPole, colour ='N Hotspots\nInduced by POLE'))+
  geom_smooth(data=df,aes(x=nmutIM341, y=nTotalHotspots - nHotspotsInducedByPole, colour = 'N Hotspots\nNot Induced by POLE'))+
  
  xlab('N Nonsynonymous Mutations in IMPACT 341 Genes')+
  ylab('N hotspots')+
  ggtitle('Observed and expected hotspot burden \nin POLE driven tumors')+
  labs(caption='plot_mut_susceptibility_info.R\nanalyze_hypothetical_and_observed_mutability_by_signature.ipynb')



ggplot(df, aes(x=nmutIM341))+
 geom_smooth(aes(y=nTotalHotspots-nHotspotsExpected, colour = 'Total N hotspots\n >than expectation'))+
 geom_smooth(aes(y=nTotalHotspots-nHotspotsInducedByPole, colour = 'N hotspots induced by POLE\n >than expectation'))+
 ylab('N hotspots > than expectation')
  
 
df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/truncatingObsExpectedPOLE.tsv',sep = '\t', header=TRUE)
#WITH PENTA CONTEXT
dfPenta <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/truncatingObsExpectedPOLEPenta.tsv',sep = '\t', header=TRUE)

ggplot()+
  geom_smooth(data = dfPenta, method='glm', aes(x=nmut, y=nOncogeneTrunc-expectedOncogeneTrunc, colour='Oncogene truncating mutations\nEXPECTATION METHOD 2\n'))+
  geom_smooth(data=df, method='glm', aes(x=nmut, y=nOncogeneTruncating-nOncogeneExpected, colour='Oncogene truncating mutations\nEXPECTATION METHOD 1\n'))+
  geom_ribbon(data=df,)
  
  geom_smooth(data = dfPenta, method='glm', aes(x=nmut, y=nTSGTrunc-expectedTSGTrunc, colour='TSG truncating mutations\nEXPECTATION METHOD 2\n'))+
  geom_smooth(data=df, method='glm', aes(x=nmut, y=nTsgTruncating-nTSGExpected, colour='TSG truncating mutations\nEXPECTATION METHOD 1\n'))+
  
  scale_color_manual(values=c('#ffb347', '#f05e23', '#00008b',  'blue'))+
  ggtitle('Truncating mutations in POLE tumors')+
  ylab('N truncating mutations > than expected')+
  xlab('N mutations in IM-341 genes')+
  labs(caption='plot_mut_susceptibility_info.R\nanalyze_hypothetical_and_observed_mutability_by_signature.ipynb')

###
#########
###############
####################
###############
#########
######
###

#MSI

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/msiObsVsExpected.tsv',sep = '\t', header=TRUE)
ggplot(df, aes(x= nIndels))+
  geom_smooth(aes(y=nTsg-nTsgExpected))+
  geom_point(aes(y=nTsg-nTsgExpected))+
  ylab('N TSG mutated above Expectation')+
  xlab('N Indels')+
  ggtitle('MSI Indels Preferentially Target TSGs')+
  labs(caption='plot_mut_susceptibility.R\nanalyze_hypothetical_and_observed_msi_mutability.ipynb')


###
#########
###############
####################
###############
#########
######
###

#TMZ

dfGlioma <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/hotspotObsExpectedTMZ.tsv',sep = '\t', header=TRUE)
ggplot(dfGlioma)+
  geom_smooth(aes(x = nmutIM341, y = nHotspotsExpected, colour='n hotspots expected'))+
  geom_smooth(aes(x = nmutIM341, y = nHotspotsExpected2, colour='n hotspots expected2'))+
  
  geom_point(aes(x = nmutIM341, y= nTotalHotspots - nHotspotsInducedBySig, colour='n hotspots\n likely not induced by TMZ'))+
  geom_point(aes(x = nmutIM341, y= nHotspotsInducedBySig, colour='n hotspots\n likely induced by TMZ'))+
  ggtitle('observed vs expected tmz glioma')


#OLD STUFF
df <- read.table('/Users/friedman/Desktop/WORK/hotspotFractionOfMutationBurden.tsv',sep = '\t', header=TRUE)

barColorPalette = c(
  "#0b6623", "#67b826",
  "#A4DBE8", "#3255A4",
  "#ffb347", '#CC5500',
  'red'
)
plottingLevelsNeutral  <- c('colorectalHyper', 'colorectalNormal', 
                            'endometrialNormal', 'endometrialHyper',
                            'gliomaNormal', 'gliomaHyper',
                            'theoreticalSusceptibility')

names(barColorPalette)= plottingLevelsNeutral

pHotspot <- ggplot(df, aes(group = simplifiedClass, x=hotspotFrac, colour=simplifiedClass))+
  scale_x_log10()+
  geom_density(size=2)

pOncogenic <- ggplot(df, aes(group = simplifiedClass, x=oncogenicFrac, colour=simplifiedClass))+
  scale_x_log10()+
  geom_density(size=2)

alignedPlot <- plot_grid(pHotspot, pOncogenic, ncol = 2)
fullP <- plot_grid(ggplot()+ggtitle('Signatures vs Observed Hotspot and Onc Mut Fractions'),
                   alignedPlot, ggplot()+labs(caption='plot_mut_susceptibility_info.R  oncogenic_mut_prob_simulation_by_gene.ipynb'),
                   nrow=3, rel_heights = c(.1,1,.1))

ggsave('~/Desktop/plot.pdf', plot=fullP,  width = 10, height = 5, units = c("in"))


#
####
########
##########
##############
#LATER version with simplified labels

emptyTheme <- theme(axis.line = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    axis.ticks.y = element_blank())

#HOTSPOTS
dfHot <- read.table('/Users/friedman/Desktop/WORK/hotspotFractionOfMutationBurden.tsv',sep = '\t', header=TRUE)


dHotspot <- dfHot[dfHot$class != 'theoreticalSusceptibility',]
maxHotspotQuadNucDf = setNames(data.frame(0.0138, 'Hotspot prob @ most hotspot philic residue: T(C>G)G'), c('hF', 'label'))
maxHotspotSignatureDf = setNames(data.frame(0.00417, 'Hotspot prob @ most hotspot philic signature: POLE'), c('hF', 'label'))
dHotspot <- bind_rows(dHotspot, maxHotspotQuadNucDf, maxHotspotSignatureDf)
pHotspot <- ggplot(dHotspot, aes(group = simplifiedClass, x=hotspotFrac, colour=simplifiedClass))+
  geom_density(size=2)+
  geom_label_repel(aes(x=hF, label=label, y=0),
                  arrow = arrow(length = unit(0.03, "npc"), type = "closed", ends = "last"),
                  size=1.5,
                  segment.alpha = 0.75,
                  segment.size  = 1,
                  nudge_y = .5,
                  force=10,
                  direction    = "y",
                  segment.color = "grey50")+
  emptyTheme+
  theme(legend.position = 'none')+
  scale_x_log10(breaks=c(0, 0.004, 0.014, 0.05, .5, 1))+
  ylab('Density')+
  xlab('Log fraction of case TMB due to hotspots')+
  ggtitle('Hotspots')+
  theme(plot.title = element_text(hjust=.5, face='bold'))


##ONCOGENIC
dfOnc <- read.table('/Users/friedman/Desktop/WORK/oncogenicFractionOfMutationBurden.tsv',sep = '\t', header=TRUE)

dOncogenic <- dfOnc[dfOnc$class != 'theoreticalSusceptibility',]
maxOncogenicQuadNucDf = setNames(data.frame(0.2, 'Oncogenic prob @ most onco-philic residue: T(C>A)G'), c('oF', 'label'))
maxOncogenicSignatureDf = setNames(data.frame(0.06, 'Oncogneic prob @ most onco-philic signature: POLE'), c('oF', 'label'))
dOncogenic <- bind_rows(dOncogenic, maxOncogenicQuadNucDf, maxOncogenicSignatureDf)
pOncogenic <- ggplot(dOncogenic, aes(group = simplifiedClass, x=oncogenicFrac, colour=simplifiedClass))+
  geom_density(size=2)+
  geom_label_repel(aes(x=oF, label=label, y=0),
                   arrow = arrow(length = unit(0.03, "npc"), type = "closed", ends = "last"),
                   size=1.5,
                   segment.alpha = 0.75,
                   segment.size  = 1,
                   nudge_y = .5,
                   force=10,
                   direction    = "y",
                   segment.color = "grey50")+
  emptyTheme+
  scale_x_log10(breaks=c(0, .01, 0.06, 0.2, .5, 1))+
  ylab('Density')+
  xlab('Log fraction of case TMB due to Oncogenic muts')+
  ggtitle('Oncogenic Muts')+
  theme(plot.title = element_text(hjust=.5, face='bold'))

alignedPlot <- plot_grid(pHotspot, pOncogenic, ncol=2, rel_widths=c(.8,1))
fullPlot <- plot_grid(ggplot()+ggtitle('Activating mutation fractions of TMB')+theme(plot.title = element_text(hjust=.5, face='bold')),
                      alignedPlot, ggplot()+labs(caption = 'plot_mut_susceptibility_info.R, oncogenic_mut_prob_simulation_by_gene.ipynb'),
                      rel_heights = c(.1,1,.05), nrow=3)

ggsave('~/Desktop/plot.pdf', plot=fullPlot,  width = 8, height = 3, units = c("in"))



#######
#
#
#
#
#
#
#################################

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/tmzAlleleCounts.tsv',sep = '\t', header=TRUE)

ggplot(df, aes(x=count))+
  geom_histogram(stat='count')+
  ggtitle('N occurences of \nTMZ inducible alleles')







  
  
  
  
  
  
  
  
  




















