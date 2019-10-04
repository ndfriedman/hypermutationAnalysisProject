#written by noah friedman
library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)
library(ggrepel)



#
#
#
####PLOT observed vs possible oncogenic mutations by genes
df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/geneObservedFracs.tsv',sep = '\t', header=TRUE)

plt <- ggplot(df, aes(x=reorder(gene, fractionOfPossibleOncMutsObserved), y=fractionOfPossibleOncMutsObserved))+
  geom_bar(stat='identity')+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=5))+
  geom_text_repel(aes(label=displayName), nudge_y=0.5)+
  ggtitle('Fraction of Possible Oncogenic SNPs observed in IMPACT')

plt2 <- ggplot(df[df$ratio < 1,], aes(x=reorder(gene, ratio), y=ratio))+
  geom_bar(stat='identity')+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=5))+
  geom_text_repel(aes(label=displayName), nudge_y=0.1)+
  ggtitle('Fraction of Possible Muts observed in IMPACT')
  
alignedPlot <- plot_grid(plt, plt2, nrow = 2)
ggsave('~/Desktop/plot.pdf', plot=alignedPlot,  width = 30, height = 10, units = c("in"))

pC <- ggplot(df, aes(x = reorder(gene, nPossibleOncogenic), y=nPossibleOncogenic))+
  geom_bar(stat='identity')+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=5))+
  geom_text_repel(aes(label=displayName), nudge_y=5)+
  ggtitle('N possible onc mutants in IMPACT')

ggsave('~/Desktop/plot.pdf', plot=pC,  width = 25, height = 5, units = c("in"))

pGC <- ggplot(df, aes(x = reorder(gene, nPossible), y=nPossible))+
  geom_bar(stat='identity')+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=5))+
  geom_text_repel(aes(label=displayName), nudge_y=125)+
  ggtitle('N possible non-synonymous mutants in IMPACT')

ggsave('~/Desktop/plot.pdf', plot=pGC,  width = 25, height = 5, units = c("in"))


#
#
#
#####PLOT observed vs possible quadnuc mutations by genes
df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/quadNucObservedFracs.tsv',sep = '\t', header=TRUE)

p1 <- ggplot(df, aes(x=reorder(quadNuc, nObservedQuadNuc/nPossibleMuts), y=nObservedQuadNuc/nPossibleMuts))+
  geom_bar(stat='identity', aes(fill=changeType))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=7))+
  scale_fill_manual(values= c('blue', 'black', 'red', 'gray', 'green', 'pink'))+
  ggtitle('Fraction of possible mutations observed')

p2 <- ggplot(df, aes(x=reorder(quadNuc, nObservedOncAtQuadNuc/nPossibleOncMuts), y=nObservedOncAtQuadNuc/nPossibleOncMuts))+
  geom_bar(stat='identity', aes(fill=changeType))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=7))+
  scale_fill_manual(values= c('blue', 'black', 'red', 'gray', 'green', 'pink'))+
  ggtitle('Fraction of possible Oncogenic Muts Observed')

caption <- ggplot() + labs(caption = 'compare_observed_and_simulated_mutation_characteristics.ipynb, plotObservedVsPossibleMutInfo.R')

alignedPlot <- plot_grid(p1, p2, caption, nrow = 3, rel_heights = c(1,1,.1))
ggsave('~/Desktop/plot.pdf', plot=alignedPlot,  width = 15, height = 10, units = c("in"))

#
#
#
#######PLOT quadnucs expected
df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/quadNucExpectedFracs.tsv',sep = '\t', header=TRUE)
plt <- ggplot(df, aes(x = reorder(quadNuc, probability), y=probability))+
  geom_bar(stat='identity', aes(fill=changeType))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=7))+
  scale_fill_manual(values= c('blue', 'black', 'red', 'gray', 'green', 'pink'))+
  ggtitle('probability of oncogenic mutation by quadnuc')+
  geom_text_repel(aes(label=motifLabel), nudge_y=0.03)

ggsave('~/Desktop/plot.pdf', plot=plt,  width = 15, height = 10, units = c("in"))


#
#
#
##################PLOT distinct mutations vs other stuff

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/mutationObservedVsPossibleMuts.tsv',sep = '\t', header=TRUE)

plt <- ggplot(df[df$nObservedDistinctNonSilentMuts/df$nPossibleNonSilentMuts < 1,], aes(x=reorder(gene, nObservedDistinctNonSilentMuts/nPossibleNonSilentMuts), y=nObservedDistinctNonSilentMuts/nPossibleNonSilentMuts))+
  geom_bar(stat='identity')+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=6))+
  geom_text_repel(aes(label=displayName), nudge_y=0.1)

ggsave('~/Desktop/plot.pdf', plot=plt,  width = 30, height = 10, units = c("in"))


#
#
#
############################PLOT distinct vs possible by variant type

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/variantTypeObsVsPossible.tsv',sep = '\t', header=TRUE)
plt <- ggplot(df, aes(x=mutationType, y=frac))+
  geom_bar(stat='identity')+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=6))+
  ggtitle('Fraction of possible mutations observed in IMPACT')

#
#
#
##############################################PLOT SIGS BY ONC MUT FRACTION

make_sigs_plot <- function(df, title, segY){
  p <- ggplot(df, aes(x=nmut, y=nOncogenic/nmut))+
    geom_point(alpha=0.15)+
    geom_segment(aes(x=0, xend=max(df$nmut), y=segY, yend=segY), colour= 'red')+
    ggtitle(title)+
    scale_x_log10()+
    geom_smooth()+
    ylim(0,1)
  return(p)
}

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/oncMutFracsBySignatures.tsv',sep = '\t', header=TRUE)
df <- df[df$nmut > 30,]
plts <- plot_grid(
  make_sigs_plot(df[df$Signature_Aetiology == 'POLE',], 'POLE', 0.065),
  make_sigs_plot(df[df$Signature_Aetiology == 'MMR',], 'MMR', 0.17),
  make_sigs_plot(df[df$Signature_Aetiology == 'TMZ',], 'TMZ', 0.045),
  make_sigs_plot(df[df$Signature_Aetiology == 'APOBEC',], 'APOBEC', 0.08),
  make_sigs_plot(df[df$Signature_Aetiology == 'UV',], 'UV', 0.05),
  make_sigs_plot(df[df$Signature_Aetiology == 'SMOKING',], 'SMOKING', 0.075),
  nrow =2, ncol=3
)

ggsave('~/Desktop/plot.pdf', plot=plts,  width = 15, height = 10, units = c("in"))

