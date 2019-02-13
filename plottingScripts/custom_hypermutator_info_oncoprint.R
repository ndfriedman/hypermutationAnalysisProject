#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

plot_oncoplot_endometrial <- function(df, tileNames){
  mutBurdenBar <- ggplot(df, aes(x=reorder(Tumor_Sample_Barcode, Nmut_Mb), y=Nmut_Mb))+
    geom_bar(stat='identity')+
    scale_y_log10()+
    ylab('Log Nmut_Mb')+
    theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.x=element_blank(),
            axis.title.x=element_blank())+
    ggtitle('Mutation landscape of Uterine Endometriod/Endometriods Tumors')+
    geom_vline(aes(xintercept=140))+ #constants based on data
    geom_vline(aes(xintercept=93))+
    geom_text(aes(x=106, y=90, label='High Mut Burden: 20-80 mut_mb'), size =6)+
    geom_text(aes(x=152, y=900, label='Ultra High Mut Burden: 80+ mut_mb'), size =6)
  ptenOncoprint <- ggplot(df, aes(x=reorder(Tumor_Sample_Barcode, Nmut_Mb), y=1, 
                                  fill=factor(PTEN_status, levels=tileNames)))+
    geom_tile()+
    get_empty_theme()+
    scale_fill_manual(values=barColorPalette, drop=FALSE)+
    theme(legend.position="none")+
    theme(axis.title.y=element_text(size=14,face="bold"))+
    ylab('PTEN')
  pik3caOncoprint <- ggplot(df, aes(x=reorder(Tumor_Sample_Barcode, Nmut_Mb), y=1, 
                                    fill=factor(PIK3CA_status, levels=tileNames)))+
    geom_tile()+
    get_empty_theme()+
    scale_fill_manual(values=barColorPalette, drop=FALSE)+
    theme(legend.position="none")+
    theme(axis.title.y=element_text(size=14,face="bold"))+
    ylab('PIK3CA')
  arid1aOncoprint <- ggplot(df, aes(x=reorder(Tumor_Sample_Barcode, Nmut_Mb), y=1, 
                                    fill=factor(ARID1A_status, levels=tileNames)))+
    geom_tile()+
    get_empty_theme()+
    scale_fill_manual(values=barColorPalette, drop=FALSE)+
    theme(legend.position="none")+
    theme(axis.title.y=element_text(size=14,face="bold"))+
    ylab('ARID1A')
  tp53Oncoprint <- ggplot(df, aes(x=reorder(Tumor_Sample_Barcode, Nmut_Mb), y=1, 
                                  fill=factor(TP53_status, levels=tileNames)))+
    geom_tile()+
    get_empty_theme()+
    scale_fill_manual(values=barColorPalette, drop=FALSE)+
    theme(legend.position="none")+
    theme(axis.title.y=element_text(size=14,face="bold"))+
    ylab('TP53')
  nf1Oncoprint <- ggplot(df, aes(x=reorder(Tumor_Sample_Barcode, Nmut_Mb), y=1, 
                                  fill=factor(NF1_status, levels=tileNames)))+
    geom_tile()+
    get_empty_theme()+
    scale_fill_manual(values=barColorPalette, drop=FALSE)+
    theme(legend.position="none")+
    theme(axis.title.y=element_text(size=14,face="bold"))+
    ylab('NF1')
  ctnnb1Oncoprint <- ggplot(df, aes(x=reorder(Tumor_Sample_Barcode, Nmut_Mb), y=1, 
                                  fill=factor(CTNNB1_status, levels=tileNames)))+
    geom_tile()+
    scale_fill_manual(values=barColorPalette, drop=FALSE)+
    get_empty_theme()+
    theme(legend.position="none")+
    theme(axis.title.y=element_text(size=14,face="bold"))+
    ylab('CTNNB1')
    
  signatureOncoprint <-ggplot(dataF, aes(x=reorder(Tumor_Sample_Barcode, Nmut_Mb), y=1, 
                                         fill=signatureAetiology))+
    geom_tile()+
    get_empty_theme()+
    theme(axis.text.x=element_text(angle = 60, hjust = .5, size=2))+theme(axis.title.y=element_text(size=14,face="bold"))+
    ylab('SIGNATURE')+
    scale_fill_manual(values=signatureColorPalette, drop=FALSE)+
    labs(caption = "data:double_mut_oncoprint.py; plotting: custom_hypermutator_info_oncoprint.R")
  
  legendSignatures = get_legend(signatureOncoprint)
  signatureOncoprint <- signatureOncoprint + theme(legend.position="none")
  
  #Throw away code that exists here only to give us the legend
  legendP <- ggplot(df, aes(x=reorder(Tumor_Sample_Barcode, Nmut_Mb), y=1, 
                                      fill=factor(CTNNB1_status, levels=tileNames)))+
                                      geom_tile()+
                                      scale_fill_manual(values=barColorPalette, drop=FALSE)
  legendP <- legendP + guides(fill=guide_legend(title="Mutation/CNA Status"))
  legend <- get_legend(legendP)
    
  #align the panels then the legend too
  alignedPlot <- plot_grid(
    mutBurdenBar, ptenOncoprint, pik3caOncoprint, arid1aOncoprint, nf1Oncoprint, tp53Oncoprint, ctnnb1Oncoprint, signatureOncoprint, 
    align='hv', nrow=8, ncol=1, rel_heights = c(1,.3,.3,.3, .3, .3, .3, .2)
  )
  legendPlot <- plot_grid(
    ggplot(), ggplot(), legend, ggplot(),ggplot(), legendSignatures, align='hv', nrow=6, ncol=1 
  )
  alignedPlotWithLegend <- plot_grid(alignedPlot, legendPlot, align='hv', nrow=1, ncol=2, rel_widths = c(1,.15))
  return(alignedPlotWithLegend)
}


barColorPalette = c("#003366", #Loh + double onco is dark blue  
                    "orange", #oncogenic plus loh orange
                    "#00008B", #double oncogenic (navy)
                    "#2A52BE", #oncogenic plus other (pale blue)
                    "red", #oncogenic only is red
                    "yellow", #loh only
                    "#D3D3D3" #neither is gray
)

signatureColorPalette = c(
  "#00FFFF", #age
  "#00DD5D", #sig14,
  "#267574", #mmr
  "white", #not enough mutations
  "gray", #other
  "#ADFF2F" #pole
)

dataF <- read.table('~/Desktop/WORK/dataForLocalPlotting/endometrialOncoprintData.tsv',sep = '\t', header=TRUE)
tNames = c('Double Oncogenic & LOH', 'Oncogenic+LOH', 'DoubleOncogenic', 'Oncogenic+Other', 'OncogenicMut', 'LohOnly', 'Neither oncogenic mut nor loh')
plt <- plot_oncoplot_endometrial(dataF, tNames)

ggsave('~/Desktop/plot.pdf', plot=plt,  width = 32, height = 20, units = c("in"))



dim(dataF[dataF$Nmut_Mb < 20,])

dataF$Tumor_Sample_Barcode
