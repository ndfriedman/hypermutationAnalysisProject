#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
require(ggpubr)
require(ggridges)
library(scales)
library(viridis)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

#SAME BOX PLOT FUNCTION BUT NOT HAMSTRUNG BY BS
vanillaBoxPlot <- function(df, xValParam, yColParam, comps=NA, compPos=NA, compYPositions=NA, maxY, title, yLabelText){
  
  df <- my.rename(df, xValParam, "xVal")
  df <- my.rename(df, yColParam, "yCol")
  
  #TODO FIX THIS AGAIN
  p<- ggplot(df, aes(reorder(xVal, -orderingVal), y=yCol)) +
    geom_boxplot(fatten = NULL)+ #TODO FIX THIS TO PLOT THE MIDDLE LINE AS THE MEAN!!!
    stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1, linetype = "solid")+
    coord_cartesian(ylim=c(0, maxY+10))+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
    ylab(yLabelText)+
    xlab(xValParam)+
    theme(legend.position="none")+
    ggtitle(title)+
    labs(caption = "data:analyze_hypermutator_hotspot_burden.py; plotting: hypermutator_driver_plottings.R")
  
  
  p<- p+ geom_jitter(shape=16, position=position_jitter(0.1))
  
  p<- p + stat_compare_means(comparisons = 
                               comps,
                             label.y=
                               compPos, tip.length=.025, label = "p.signif")
  return(p)
}

plotBoxPlot <- function(df, xValParam, yColParam, comps=NA, compPos=NA, compYPositions=NA, maxY, title, yLabelText){
  
  df <- my.rename(df, xValParam, "xVal")
  df <- my.rename(df, yColParam, "yCol")

  #TODO FIX THIS AGAIN
  p<- ggplot(df, aes(reorder(xVal, -orderingVal), y=yCol)) +
    geom_boxplot(fatten = NULL)+ #TODO FIX THIS TO PLOT THE MIDDLE LINE AS THE MEAN!!!
    stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1, linetype = "solid")+
    coord_cartesian(ylim=c(0, maxY+10))+
    #theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
    theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(), axis.title.x=element_blank())+
    ylab(yLabelText)+
    
    #ADD EXTRA TEXT LABELS FOR READABILITY
    #OLD TEXT LOCATIONS
    #geom_text(aes(x=2, y=-2, label='Non-small cell lung cancer', size=20), fontface="italic")+
    #geom_text(aes(x=5, y=-2, label='Bladder cancer', size=20), fontface="italic")+
    #geom_text(aes(x=8, y=-2, label='Melanoma', size=20), fontface="italic")+
    #geom_text(aes(x=11, y=-2, label='Glioma', size=20), fontface="italic")+
    #geom_text(aes(x=14, y=-2, label='Colorectal cancer', size=20), fontface="italic")+
    #geom_text(aes(x=17, y=-2, label='Endometrial cancer', size=20), fontface="italic")+
    
    #NEW TEXT LOCATIONS
    geom_text(aes(x=1.5, y=-2, label='Non-small cell lung cancer', size=20), fontface="italic")+
    geom_text(aes(x=3.5, y=-2, label='Bladder cancer', size=20), fontface="italic")+
    geom_text(aes(x=5.5, y=-2, label='Melanoma', size=20), fontface="italic")+
    geom_text(aes(x=7.5, y=-2, label='Glioma', size=20), fontface="italic")+
    geom_text(aes(x=9.5, y=-2, label='Colorectal cancer', size=20), fontface="italic")+
    geom_text(aes(x=12, y=-2, label='Endometrial cancer', size=20), fontface="italic")+
    geom_text(aes(x=14.5, y=-2, label='Other Hypermutant', size=20), fontface="italic")+
    
    theme(legend.position="none")
  
  #p<- p + stat_compare_means(comparisons = comps, label.y=compYPositions, tip.length=.01)#
  p<- p+ geom_jitter(shape=16, position=position_jitter(0.1))
  
  p<- p + stat_compare_means(comparisons = 
                               comps,
                             label.y=
                               compPos, tip.length=.025, label = "p.signif")
  
  p <- p + ggtitle('Landscape of Oncogenic burden across hypermutated and non-hypermutated cancers')
  
  return(p)
}


plot_box_plt_with_panels <- function(df, Comps, CompPos){
  boxPlt <- plotBoxPlot(df, 'cohort', 'nOncogenicMutations', 
                        maxY=max(df$nOncogenicMutations), title='N Oncogenic Mutations', yLabelText = 'n oncogenic mutations',
                        comps=Comps, compPos=CompPos)
  
  nCasesPanel <- ggplot(df, aes(reorder(cohort, -orderingVal), fill=cancer_type))+
    geom_bar(stat='count')+
    scale_y_log10()+
    theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(), axis.title.x=element_blank(), axis.line = element_blank())+
    theme(legend.position="none")+
    theme(axis.title.y=element_text(angle = 45, hjust = .2, size=14, face="bold"))+
    ylab('N cases in cohort')
  
  averageNMutPanel <- ggplot(df, aes(x=reorder(cohort, -orderingVal), y=Nmut_Mb))+
    stat_summary(fun.y = mean, geom = "point", aes(ymax = ..y.., ymin = ..y..),
    size = 6, colour='Red')+
    theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(), axis.title.x=element_blank(), axis.line = element_blank())+
    theme(axis.title.y=element_text(angle = 45, hjust = .2, size=14, face="bold"))+
    ylab('Average mut_mb in cohort')
  
  signatureOncoprint <-ggplot(df, aes(x=reorder(cohort, -orderingVal), 
                                         fill=signatureAetiology))+
    geom_bar(position="fill")+
    get_empty_theme()+
    theme(axis.title.y=element_text(angle = 45, hjust = .2, size=14, face="bold"))+
    ylab('SIGNATURE')+
    theme(axis.text.x=element_text(angle = 90, hjust = .5, size=11, face="bold"))+theme(axis.title.y=element_text(size=14,face="bold"))+
    scale_fill_manual(values=signatureColorPalette, drop=FALSE)+
    labs(caption = "data:analyze_hypermutator_hotspot_burden.py; plotting: hypermutator_driver_plottings.R")
  
  signatureOncoprint <- signatureOncoprint + theme(legend.position="none")
  legendSigs <- get_legend(signatureOncoprint)

  #ALIGN THE PLOTS
  alignedPlot <- plot_grid(
    boxPlt, nCasesPanel, averageNMutPanel, signatureOncoprint,
    align='v',nrow=4, ncol=1, rel_heights = c(.6,.1, .075, .2)
  )
  
  #ALIGN THE LEGENDS
  legendPlot <- plot_grid(
    ggplot(),
    legendSigs,
    ggplot(),
    align = 'hv', nrow = 2, ncol= 1, rel_heights = c(1,.2, .1)
  )
  #ALIGN EVERYTHING
  alignedPlotWithLegend <-
    plot_grid(
      alignedPlot, legendPlot,
      align = 'h', nrow = 1, ncol= 2, rel_widths = c(1,.2)
    )
  
  return(alignedPlotWithLegend)
}

plot_signature_bar <- function(df, cancerType){
  plt <- ggplot(df[df$cancer_type == cancerType,],
                aes(x=reorder(Tumor_Sample_Barcode, Nmut), y=1, 
                    fill=signatureAetiology))+
    geom_tile()+
    get_empty_theme()+
    ylab('SIGNATURE')+
    scale_fill_manual(values=signatureColorPalette, drop=FALSE)+
    get_empty_theme()+
    ylab(cancerType)+
    theme(axis.title.y=element_text(angle = 0, hjust = .2, size=14, face="bold"))
}

plot_signature_information <- function(df){
  alginedPlot <- plot_grid(
    plot_signature_bar(df, 'Non-Small Cell Lung Cancer')+theme(legend.position = "none"),
    plot_signature_bar(df, 'Bladder Cancer')+theme(legend.position = "none"),
    plot_signature_bar(df, 'Melanoma')+theme(legend.position = "none"),
    plot_signature_bar(df, 'Glioma')+theme(legend.position = "none"),
    plot_signature_bar(df, 'Colorectal Cancer')+theme(legend.position = "none"),
    plot_signature_bar(df, 'Endometrial Cancer')+theme(legend.position = "none"),
    nrow=6, align='v'
  )
  legend <- get_legend(plot_signature_bar(df, 'Glioma'))
  plotWithLegend <- plot_grid(alginedPlot, legend, ncol=2, rel_widths = c(1,.2))
  return(plotWithLegend)
}

plot_signature_mut_burden_summary <- function(df){
  plt <- ggplot(df, aes(x=reorder(dominantSignature, -orderingVal), y=Nmut_Mb))+
    geom_boxplot(fatten = NULL, alpha=0.15)+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
    stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1, linetype = "solid")+
    scale_y_continuous(breaks = c(0,10,20, 50,100,500,1000,2000))+
    coord_trans(y='log10')+
    geom_jitter(shape=16, position=position_jitter(0.1), alpha=0.15)+
    ggtitle('Mutation burden distribution by signatures in IMPACT (cases where we can infer signatures only)')+
    xlab('Signature')+
    labs(caption = "data:double_mutation_vaf_analysis.py; plotting: hypermutator_driver_plotting.R")
  return(plt)
}

#THE NEWEST VERSION OF THIS PLOTTING TOOL WHICH MAKES box plots where we have defined the signature cohorts in advance
plot_box_plt_with_panels_sig_cohorts_defined <- function(df, Comps, CompPos){
  boxPlt <- plotBoxPlot(df, 'cohort', 'nOncogenicMutations', 
                        maxY=max(df$nOncogenicMutations), title='N Oncogenic Mutations', yLabelText = 'n oncogenic mutations',
                        comps=Comps, compPos=CompPos)
  
  nCasesPanel <- ggplot(df, aes(reorder(cohort, -orderingVal), fill=cancer_type_fill))+
    geom_bar(stat='count')+
    scale_y_log10()+
    theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(), axis.title.x=element_blank(), axis.line = element_blank())+
    theme(legend.position="none")+
    theme(axis.title.y=element_text(angle = 45, hjust = .2, size=14, face="bold"))+
    ylab('N cases in cohort')
  
  averageNMutPanel <- ggplot(df, aes(x=reorder(cohort, -orderingVal), y=Nmut_Mb))+
    stat_summary(fun.y = mean, geom = "bar", aes(ymax = ..y.., ymin = ..y..),
                 size = 6, fill='orange')+
    stat_summary(fun.data=mean_se, geom = "errorbar", width = 0.25, 
                 size = 2, col = "red")+
    theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(), axis.title.x=element_blank(), axis.line = element_blank())+
    ylab('Average mut_mb in cohort')
  
  averageAgePanel <- ggplot(df, aes(reorder(cohort, -orderingVal), y=ageAtSequencing))+
    stat_summary(fun.data=mean_se, geom = "errorbar", width = 0.25, 
                 size = 2, col = "grey57") +
    stat_summary(fun.y = mean, geom = "point", aes(ymax = ..y.., ymin = ..y..),
                 size = 3, colour='Black')+
    theme(axis.title.y=element_text(angle = 45, hjust = .2, size=14, face="bold"))+
    theme(axis.text.x=element_text(angle = 90, vjust = .01, size=15, face="bold"))+
    xlab('Cohort')+
    theme(legend.position="none")+
    labs(caption = "data:analyze_hypermutator_hotspot_burden.py; plotting: hypermutator_driver_plottings.R")
  
  
  
  #ALIGN THE PLOTS
  alignedPlot <- plot_grid(
    boxPlt, nCasesPanel, averageNMutPanel, averageAgePanel, 
    align='v',nrow=4, ncol=1, rel_heights = c(.5,.1, .1, .25)
  )
  
  return(alignedPlot)
}



#CODE TO PLOT EVERYTHING WITH CURVES
#signatureAetiology
plot_density_table_panel <- function(df, cancerType){
  density <- ggplot(df[df$cancer_type == cancerType,], aes(x=Nmut_Mb, y=cancer_type))+
    geom_density_ridges()+
    scale_x_continuous(trans = log2_trans(), breaks=c(5,10, 20, 50,100,500))+
    theme(axis.text.x = element_blank(), axis.ticks.x  = element_blank(),
    axis.title.x=element_blank(),axis.line = element_blank())+
    ylab('')
  pltJitter <- ggplot(df1[df1$cancer_type == cancerType,], aes(x=Nmut_Mb, y=cancerType, colour=signatureAetiology)) +
    geom_jitter(shape=16, position=position_jitter(0.1))+
    scale_colour_manual(values=signatureColorPalette, drop=FALSE)+
    scale_x_continuous(trans = log2_trans(), breaks=c(5,10, 20, 50,100,500))+
    theme(legend.position="none", axis.line = element_blank())+
    theme(plot.margin=unit(c(-5,1,1,1), "in"))+
    ylab('')
  alignedPlot <- plot_grid(density, pltJitter, nrow=2, rel_heights = c(1,.1))
  return(alignedPlot)
}

#there is some difficulty in properly plotting this
plot_full_density_tables <- function(df){
  return(plot_grid(
    plot_density_table_panel(df, 'Melanoma'),
    plot_density_table_panel(df, 'Glioma'),
    plot_density_table_panel(df, 'Colorectal Cancer'),
    plot_density_table_panel(df, 'Endometrial Cancer'),
    nrow=2, ncol=2
  ))
}

plot_density_tables_sans_signatures <- function(df, cancerType){
  densityPlot <- ggplot(df[df$cancer_type%in%c('Endometrial Cancer', 'Colorectal Cancer', 'Glioma', 'Breast Cancer', 'Melanoma', 'Bladder Cancer', 'Prostate Cancer', 'Non-Small Cell Lung Cancer', 'Pancreatic Cancer', 'Soft Tissue Sarcoma'),], aes(x=Nmut_Mb, y=cancer_type))+
    geom_density_ridges(alpha=0.9)+
    scale_x_continuous(trans = log2_trans(), breaks=c(5,10, 20, 50,100,500))+
    ggtitle('Mutation burden distribution by cancer type')+
    labs(caption = "data:analyze_hypermutator_hotspot_burden.py; plotting: hypermutator_driver_plottings.R")
  return(densityPlot)
}

df1 <- read.table('~/Desktop/WORK/dataForLocalPlotting/signaturePlotting.tsv', sep = '\t', header=TRUE)
#plt <- plot_full_density_tables(df1)
plt <- plot_density_tables_sans_signatures(df1)

ggsave('~/Desktop/plt.pdf', plot=plt, width = 20, height = 20, units = c("in"), limitsize = FALSE)

############################NEW VERSION

df <- read.table('~/Desktop/WORK/dataForLocalPlotting/mutburdenBoxplotV2.tsv', sep = '\t', header=TRUE)

#comparissons <- list( 
#  c("Non-Small Cell Lung Cancer__not_high", 'Non-Small Cell Lung Cancer__high'), c("Non-Small Cell Lung Cancer__high", "Non-Small Cell Lung Cancer__hypermutant"),c("Non-Small Cell Lung Cancer__not_high", "Non-Small Cell Lung Cancer__hypermutant"),
#  c("Bladder Cancer__not_high", 'Bladder Cancer__high'), c("Bladder Cancer__high", "Bladder Cancer__hypermutant"),c("Bladder Cancer__not_high", "Bladder Cancer__hypermutant"),
#  c("Melanoma__not_high", 'Melanoma__high'), c("Melanoma__high", "Melanoma__hypermutant"),c("Melanoma__not_high", "Melanoma__hypermutant"),
#  c("Glioma__not_high", 'Glioma__high'), c("Glioma__high", "Glioma__hypermutant"),c("Glioma__not_high", "Glioma__hypermutant"),
#  c("Colorectal Cancer__not_high", 'Colorectal Cancer__high'), c("Colorectal Cancer__high", "Colorectal Cancer__hypermutant"),c("Colorectal Cancer__not_high", "Colorectal Cancer__hypermutant"),
#  c("Endometrial Cancer__not_high", 'Endometrial Cancer__high'), c("Endometrial Cancer__high", "Endometrial Cancer__hypermutant"),c("Endometrial Cancer__not_high", "Endometrial Cancer__hypermutant"))

#comparissonPositions <- c(25,30,35,
#                          25,60,65,
#                          25,35,40,
#                          25,65,70,
#                          35,70,75,
#                          40,80,85)
comparissons<- NA
comparissonPositions <-NA

signatureColorPalette = c(
  "#00FFFF", #age
  "red", #APOBEC
  "#FF1493", #BRCA
  "#00DD5D", #sig14,
  "#267574", #mmr
  "#F2F2F2", #not enough mutations
  "pink", #other
  "#ADFF2F", #pole
  "orange", #smoking
  "dark blue", #TMZ
  "yellow" #UV
)

#plt <- plot_box_plt_with_panels(df, comparissons, comparissonPositions)
plt <- plot_box_plt_with_panels_sig_cohorts_defined(df, comparissons, comparissonPositions)
ggsave('~/Desktop/plt.pdf', plot=plt, width = 30, height = 25, units = c("in"), limitsize = FALSE)





#$#$#$#NOAH 
#$#$#$#$NOAH 
#$$#$#$#NOAH
df = read.table('~/Desktop/WORK/dataForLocalPlotting/vanillaImpactSigsFile.tsv', sep = '\t', header=TRUE)
plot_density_table(df1)


#OTHER PLOTS BESIDES THE BIG ONE

comparissonPositionsAge <- c(90,95,100,
                             90,95,100,
                             90,95,100,
                             90,95,100,
                             90,95,100,
                             90,95,100)

pltWithAge <- vanillaBoxPlot(df, 'cohort', 'ageAtSequencing', 
            maxY=105, title='Age at Sequencing', yLabelText = 'Age at Sequencing',
            comps=comparissons, compPos=comparissonPositionsAge)

ggsave('~/Desktop/plt.pdf', plot=pltWithAge, width = 10, height = 8, units = c("in"), limitsize = FALSE)


comparissonPositionsNDoublets <- c(5,7,9,
                                   5,7,9,
                                   5,7,9,
                                   5,7,9,
                                   7,9,12,
                                   7,9,15)

pltWithNDoublets <- vanillaBoxPlot(df, 'cohort', 'nGenesWithDoubleOncogenic', 
                                   maxY=10, title='N double oncogenic mutated genes in case', yLabelText = 'N double oncogenic mutations',
                                   comps=comparissons, compPos=comparissonPositionsNDoublets)

ggsave('~/Desktop/plt.pdf', plot=pltWithNDoublets, width = 10, height = 8, units = c("in"), limitsize = FALSE)



################MAKE A PLOT THAT COMPARES SIGNATURES PREVALENCES ACROSS MUT BURDEN BANDS
df <- read.table('~/Desktop/WORK/dataForLocalPlotting/signaturePlotting.tsv', sep = '\t', header=TRUE)
#df <- df[df$Nmut > 25,]
plt <- plot_signature_information(df)

ggsave('~/Desktop/plot.pdf', plot=plt,  width = 20, height = 15, units = c("in"))




####################################
df <- read.table('~/Desktop/WORK/dataForLocalPlotting/signatureSummary.tsv', sep = '\t', header=TRUE)
plt <- plot_signature_mut_burden_summary(df)
ggsave('~/Desktop/plot.pdf', plot=plt,  width = 10, height = 7, units = c("in"))


