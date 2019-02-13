#written by Noah Friedman (a template for scripts to be excuted in the spyder environment
import sys
import argparse
import os
import pandas as pd
import numpy as np

from collections import Counter

pathPrefix = ''
if os.getcwd() == '/Users/friedman/Desktop/mnt':
	pathPrefix = '/Users/friedman/Desktop/mnt'

sys.path.append(pathPrefix + '/ifs/work/taylorlab/friedman/myUtils')
import analysis_utils 
import mutationSigUtils 
import maf_analysis_utils
import double_mutation_analysis_util

import scipy.stats

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h


#returns proto dataframe information with averages of each column
def process_mut_summary_maf_for_averages(summaryMaf, genes, suffix='_doubleOncogenic', labelMode='hypermutator'):
    
    listOfDicts = [] #the list of dictionaries we will return from this function
    sortedArr, allValues, avg = summarize_top_double_mutated_genes(summaryMaf, suffix)
    m, lowerBound, upperBound = mean_confidence_interval(allValues, confidence=0.95)
    listOfDicts.append({'label': labelMode + '_all_genes_average', 'value': avg, 'lowerConf':lowerBound, 'upperConf':upperBound})
    for gene in genes:
        col = gene + suffix
        avg = np.nanmean(summaryMaf[col])
        m, lowerBound, upperBound = mean_confidence_interval(list(summaryMaf[col]), confidence=0.95)
        listOfDicts.append({'label': labelMode + '_' + gene, 'value': avg, 'lowerConf':lowerBound, 'upperConf':upperBound})
    return listOfDicts

#a utility function that summarizes the genes that are most recurrenctly double mutated across a cohort
#returns a sorted list of the genes, all 1s/0s for what happens for each level, the average
def summarize_top_double_mutated_genes(sMaf, suffix='_doubleOncogenic'):
    tupleArr = []
    valueArr = []
    arr = []
    for column in sMaf.columns.values:
        if suffix in column:
            arr.append(column)
    for val in arr: 
        tupleArr.append((val, np.nanmean(sMaf[val])))
        valueArr = valueArr + list(sMaf[val])
    return sorted(tupleArr, key= lambda x: x[1], reverse=True), valueArr, np.nanmean([i[1] for i in tupleArr])

def assign_ordering_val(row, orderingDict):
    geneName = row['label'].split('_')[len(row['label'].split('_')) - 1] #get the substring after the last underscore (the gene name or the word "average")
    isHypermutator = True
    if row['label'][0] == 'N':
        isHypermutator = False  #the string starts with NOT-hypermutator
    if geneName == 'average':
        if isHypermutator:
            return 1000 #big value to make sure these stay at the end
        else:
            return 1001
    else:
        if isHypermutator:
            return orderingDict[geneName]
        else:
            return orderingDict[geneName] + .1

#big function that takes a cohort maf and generates statistics on doublet prevalence of nDoublesToPlot doublets     
def create_hypermutation_prevalence_summary(mutSummaryMaf, nDoublesToPlot = 10, doubleSuffix='_doubleOncogenic'):
    mutSummaryMafHypermutators = mutSummaryMaf[mutSummaryMaf['isHypermutator'] == True]
    mutSummaryMafNonHypermutators = mutSummaryMaf[mutSummaryMaf['isHypermutator'] == False]
    
    topHypermutatorGenesArr, v1, v2 = summarize_top_double_mutated_genes(mutSummaryMafHypermutators, suffix=doubleSuffix)
    topNDoublesDict = dict([(topHypermutatorGenesArr[:nDoublesToPlot][i][0].split('_')[0], i) for i in range(nDoublesToPlot)]) #complicated line that maps the first substring of the tuple to their index in a sorted list
    genesToPlot = topNDoublesDict.keys()
    
    summaryInfoHypermutatorMaf = process_mut_summary_maf_for_averages(mutSummaryMafHypermutators, genesToPlot, labelMode='hypermutator', suffix=doubleSuffix)
    summaryInfoNonHypermutatorMaf = process_mut_summary_maf_for_averages(mutSummaryMafNonHypermutators, genesToPlot, labelMode='NOT-hypermutator', suffix=doubleSuffix)
    allListOfDicts = summaryInfoHypermutatorMaf + summaryInfoNonHypermutatorMaf
    df = pd.DataFrame(allListOfDicts)
    df['orderingVal'] = df.apply(lambda row: assign_ordering_val(row, topNDoublesDict), axis=1)
    df['isHypermutation'] = df['label'].apply(lambda x: False if 'NOT-' in x else True)
    
    return df




maf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/hypermutationAnalysisProjectData/hypermutationAnalysisProjectMaf.tsv')
sigs = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/hypermutationAnalysisProjectData/hypermutationAnalysisProjectSignatures.tsv')

sigs['pid'] = sigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
sigs['cancer_type'] = sigs['pid'].apply(lambda x: cDict[x] if x in cDict else None)

hypermuationThreshold = 80

hypermutatorIds = set(sigs[sigs['Nmut_Mb'] > hypermuationThreshold]['Tumor_Sample_Barcode'])
maf['isHypermutator'] = maf['Tumor_Sample_Barcode'].apply(lambda x: True if x in hypermutatorIds else False)

#exclude intronic variants
mafSansIntrons = maf[maf['Consequence'] != 'intron_variant']

#########ENDOMETRIAL CANCER
endometrialMaf = mafSansIntrons[mafSansIntrons['cancer_type'] == 'Endometrial Cancer']
mutSummaryMafEndometrial = double_mutation_analysis_util.create_double_mutation_summary_vanilla_maf(endometrialMaf)

endometrialDoubleSummary = create_hypermutation_prevalence_summary(mutSummaryMafEndometrial, nDoublesToPlot = 10, doubleSuffix='_double')
endometrialDoubleWithOncogenicSummary =create_hypermutation_prevalence_summary(mutSummaryMafEndometrial, nDoublesToPlot = 10, doubleSuffix='_doubleWithOncogenic')
endometrialDoubleOncogenicSummary = create_hypermutation_prevalence_summary(mutSummaryMafEndometrial, nDoublesToPlot = 10, doubleSuffix='_doubleOncogenic')

endometrialDoubleSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/endometrial_double_mutation_incidence.tsv', sep='\t', index=False)
endometrialDoubleWithOncogenicSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/endometrial_double_with_oncogenic_mutation_incidence.tsv', sep='\t', index=False)
endometrialDoubleOncogenicSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/endometrial_double_with_both_oncogenic_mutation_incidence.tsv', sep='\t', index=False)

##########COLORECTAL CANCER
colorectalMaf = mafSansIntrons[mafSansIntrons['cancer_type'] == 'Colorectal Cancer']
mutSummaryMafColorectal = double_mutation_analysis_util.create_double_mutation_summary_vanilla_maf(colorectalMaf)

colorectalDoubleSummary = create_hypermutation_prevalence_summary(mutSummaryMafColorectal, nDoublesToPlot = 10, doubleSuffix='_double')
colorectalDoubleWithOncogenicSummary =create_hypermutation_prevalence_summary(mutSummaryMafColorectal, nDoublesToPlot = 10, doubleSuffix='_doubleWithOncogenic')
colorectalDoubleOncogenicSummary = create_hypermutation_prevalence_summary(mutSummaryMafColorectal, nDoublesToPlot = 10, doubleSuffix='_doubleOncogenic')

colorectalDoubleSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/colorectal_double_mutation_incidence.tsv', sep='\t', index=False)
colorectalDoubleWithOncogenicSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/colorectal_double_with_oncogenic_mutation_incidence.tsv', sep='\t', index=False)
colorectalDoubleOncogenicSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/colorectal_double_with_both_oncogenic_mutation_incidence.tsv', sep='\t', index=False)

##########GLIOMA
gliomaMaf = mafSansIntrons[mafSansIntrons['cancer_type'] == 'Glioma']
mutSummaryMafGlioma = double_mutation_analysis_util.create_double_mutation_summary_vanilla_maf(gliomaMaf)

gliomaDoubleSummary = create_hypermutation_prevalence_summary(mutSummaryMafGlioma, nDoublesToPlot = 10, doubleSuffix='_double')
gliomaDoubleWithOncogenicSummary =create_hypermutation_prevalence_summary(mutSummaryMafGlioma, nDoublesToPlot = 10, doubleSuffix='_doubleWithOncogenic')
gliomaDoubleOncogenicSummary = create_hypermutation_prevalence_summary(mutSummaryMafGlioma, nDoublesToPlot = 10, doubleSuffix='_doubleOncogenic')

gliomaDoubleSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/glioma_double_mutation_incidence.tsv', sep='\t', index=False)
gliomaDoubleWithOncogenicSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/glioma_double_with_oncogenic_mutation_incidence.tsv', sep='\t', index=False)
gliomaDoubleOncogenicSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/glioma_double_with_both_oncogenic_mutation_incidence.tsv', sep='\t', index=False)

###############BLADDER CANCER

bladderMaf = mafSansIntrons[mafSansIntrons['cancer_type'] == 'Bladder Cancer']
mutSummaryMafBladder = double_mutation_analysis_util.create_double_mutation_summary_vanilla_maf(bladderMaf)

bladderDoubleSummary = create_hypermutation_prevalence_summary(mutSummaryMafBladder, nDoublesToPlot = 10, doubleSuffix='_double')
bladderDoubleWithOncogenicSummary =create_hypermutation_prevalence_summary(mutSummaryMafBladder, nDoublesToPlot = 10, doubleSuffix='_doubleWithOncogenic')
bladderDoubleOncogenicSummary = create_hypermutation_prevalence_summary(mutSummaryMafBladder, nDoublesToPlot = 10, doubleSuffix='_doubleOncogenic')

bladderDoubleSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/bladder_double_mutation_incidence.tsv', sep='\t', index=False)
bladderDoubleWithOncogenicSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/bladder_double_with_oncogenic_mutation_incidence.tsv', sep='\t', index=False)
bladderDoubleOncogenicSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/bladder_double_with_both_oncogenic_mutation_incidence.tsv', sep='\t', index=False)

################MELANOMA

melanomaMaf = mafSansIntrons[mafSansIntrons['cancer_type'] == 'Melanoma']
mutSummaryMafMelanoma = double_mutation_analysis_util.create_double_mutation_summary_vanilla_maf(melanomaMaf)

melanomaDoubleSummary = create_hypermutation_prevalence_summary(mutSummaryMafMelanoma, nDoublesToPlot = 10, doubleSuffix='_double')
melanomaDoubleWithOncogenicSummary =create_hypermutation_prevalence_summary(mutSummaryMafMelanoma, nDoublesToPlot = 10, doubleSuffix='_doubleWithOncogenic')
melanomaDoubleOncogenicSummary = create_hypermutation_prevalence_summary(mutSummaryMafMelanoma, nDoublesToPlot = 10, doubleSuffix='_doubleOncogenic')

melanomaDoubleSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/melanoma_double_mutation_incidence.tsv', sep='\t', index=False)
melanomaDoubleWithOncogenicSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/melanoma_double_with_oncogenic_mutation_incidence.tsv', sep='\t', index=False)
melanomaDoubleOncogenicSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/melanoma_double_with_both_oncogenic_mutation_incidence.tsv', sep='\t', index=False)





print set(['TCAT', 'TTAT', 'GCAT', 'CCAT', 'CTAT', 'ATAT', 'ACAT', 'ACTA', 'CCTA', 'GCTA', 'TCTA',
           'TCAA', 'TTAA', 'ACTT', 'TCAG', 'TTAG', 'ACTC', 'ACAA', 'ATAA', 'ACAC', 'ATAC',
           'ACAG', 'ATAG', 'ACAT', 'ATAT', 'ACTT', 'TCTG', 'TCTA', 'ACGA', 'ATGA', 'ACGC', 'ATGC',
           'ACGG', 'ATGG', 'ACGT', 'ATGT', 'TTCT', 'TTCG', 'TTCC', 'TTCA'])



