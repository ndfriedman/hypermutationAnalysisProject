#written by Noah Friedman (a template for scripts to be excuted in the spyder environment
#a script used to compare the prevalence of mutations between hypermutators and non hypermutators
import sys
import argparse
import os
import pandas as pd
import numpy as np
import math

from collections import Counter

pathPrefix = ''
if os.getcwd() == '/Users/friedman/Desktop/mnt':
	pathPrefix = '/Users/friedman/Desktop/mnt'

sys.path.append(pathPrefix + '/ifs/work/taylorlab/friedman/myUtils')
import analysis_utils 
import mutationSigUtils 
import maf_analysis_utils

def make_gene_based_df(hypermutatorGeneFracDict, hypermutatorOccurenceCounter, nonHypermutatorGeneFracDict, nonHypermutatorOccurenceCounter, rankingDict, geneLengthDict):
    listOfDicts = []
    cntr = 0
    hypermutatorOccurenceDict = dict(hypermutatorOccurenceCounter)
    nonHypermutatorOccurenceDict = dict(nonHypermutatorOccurenceCounter)
    genesToIterateThrough = dict(Counter(dict(nonHypermutatorOccurenceCounter)).most_common(50)).keys()
    for gene in genesToIterateThrough:
        cntr += 1
        localD = {}
        localD['gene'] = gene 
        if gene in hypermutatorGeneFracDict:
            localD['hypermutatorGeneFrac'] = hypermutatorGeneFracDict[gene] 
        else:
            localD['hypermutatorGeneFrac'] = 0
        if gene in hypermutatorOccurenceDict:
            localD['hypermutatorOccurence'] = hypermutatorOccurenceDict[gene] 
        else:
            localD['hypermutatorOccurence'] = 0
        if gene in nonHypermutatorGeneFracDict:
            localD['nonHypermutatorGeneFrac'] = nonHypermutatorGeneFracDict[gene] 
        else:
            localD['nonHypermutatorGeneFrac'] = 0 
        if gene in nonHypermutatorOccurenceDict:
            localD['nonHypermutatorOccurence'] = nonHypermutatorOccurenceDict[gene] 
        else:
            localD['nonHypermutatorOccurence'] = 0 
        if gene in rankingDict:
            localD['orderingVal'] = rankingDict[gene]
        else:
            localD['orderingVal'] = 500
        if gene in geneLengthDict :
            localD['LogGeneLength'] = math.log(1.0*geneLengthDict[gene], 2)
        else:
            localD['LogGeneLength'] = None
        listOfDicts.append(localD)
    return pd.DataFrame(listOfDicts)
        


maf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/annotatedOncoPlusHotspotMafAllImpact_trinuc')
mafAdj = maf_analysis_utils.add_mut_effect_summary_col(maf)
mafAdj['pid'] = mafAdj['Tumor_Sample_Barcode'].apply(lambda x: x[:9])

impactSigs = pd.read_table(pathPrefix + '/ifs/res/taylorlab/impact_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
impactSigs['pid'] = impactSigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
impactSigs['cancer_type'] = impactSigs['pid'].apply(lambda x: cDict[x] if x in cDict else None)

cdsGeneLengthDict = analysis_utils.get_cds_size_targeted_by_impact(infoFilePath = pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/impact_gene_reference_signatures.tsv')

hypermutationNmutMbThresh = 80

#COLORECTAL 
colorectalSigs = impactSigs[impactSigs['cancer_type'] == 'Colorectal Cancer']
colorectalHypermutationIds = set(colorectalSigs[colorectalSigs['Nmut_Mb'] > hypermutationNmutMbThresh]['Tumor_Sample_Barcode'])
colorectalNotHypermutationIds = set(colorectalSigs[colorectalSigs['Nmut_Mb'] <= hypermutationNmutMbThresh]['Tumor_Sample_Barcode'])
colorectalHypermutatorMaf = mafAdj[mafAdj['Tumor_Sample_Barcode'].isin(colorectalHypermutationIds)]
colorectalNONHypermutatorMaf = mafAdj[mafAdj['Tumor_Sample_Barcode'].isin(colorectalNotHypermutationIds)]
occurenceCounterHypermutators, rankingDictHypermutators, fractionalDictHypermutators = maf_analysis_utils.enumerate_top_n_oncogenic_mutated_genes_across_cohort(colorectalHypermutatorMaf, n=None)
occurenceCounterNONHypermutators, rankingDictNONHypermutators, fractionalDictNONHypermutators = maf_analysis_utils.enumerate_top_n_oncogenic_mutated_genes_across_cohort(colorectalNONHypermutatorMaf, n=None)

colorectalImportantGenes = {'APC', 'KRAS', 'PIK3CA', 'SMAD4', 'TP53', 'BRAF'}
colorectalGeneDf = make_gene_based_df(fractionalDictHypermutators, occurenceCounterHypermutators, fractionalDictNONHypermutators, occurenceCounterNONHypermutators, rankingDictHypermutators, cdsGeneLengthDict)
colorectalGeneDf.to_csv('~/Desktop/WORK/dataForLocalPlotting/colorectalHypermutationComparisson.tsv', sep='\t', index=False)

#ENDOMETRIAL
endometrialSigs = impactSigs[impactSigs['cancer_type'] == 'Endometrial Cancer']
endometrialHypermutationIds = set(endometrialSigs[endometrialSigs['Nmut_Mb'] > hypermutationNmutMbThresh]['Tumor_Sample_Barcode'])
endometrialNotHypermutationIds = set(endometrialSigs[endometrialSigs['Nmut_Mb'] <= hypermutationNmutMbThresh]['Tumor_Sample_Barcode'])
endometrialHypermutatorMaf = mafAdj[mafAdj['Tumor_Sample_Barcode'].isin(endometrialHypermutationIds)]
endometrialNONHypermutatorMaf = mafAdj[mafAdj['Tumor_Sample_Barcode'].isin(endometrialNotHypermutationIds)]
occurenceCounterHypermutators, rankingDictHypermutators, fractionalDictHypermutators = maf_analysis_utils.enumerate_top_n_oncogenic_mutated_genes_across_cohort(endometrialHypermutatorMaf, n=None)
occurenceCounterNONHypermutators, rankingDictNONHypermutators, fractionalDictNONHypermutators = maf_analysis_utils.enumerate_top_n_oncogenic_mutated_genes_across_cohort(endometrialNONHypermutatorMaf, n=None)

endometrialGeneDf = make_gene_based_df(fractionalDictHypermutators, occurenceCounterHypermutators, fractionalDictNONHypermutators, occurenceCounterNONHypermutators, rankingDictHypermutators, cdsGeneLengthDict)
endometrialGeneDf.to_csv('~/Desktop/WORK/dataForLocalPlotting/endometrialHypermutationComparisson.tsv', sep='\t', index=False)

#Glioma
gliomaSigs = impactSigs[impactSigs['cancer_type'] == 'Glioma']
gliomaHypermutationIds = set(gliomaSigs[gliomaSigs['Nmut_Mb'] > hypermutationNmutMbThresh]['Tumor_Sample_Barcode'])
gliomaNotHypermutationIds = set(gliomaSigs[gliomaSigs['Nmut_Mb'] <= hypermutationNmutMbThresh]['Tumor_Sample_Barcode'])
gliomaHypermutatorMaf = mafAdj[mafAdj['Tumor_Sample_Barcode'].isin(gliomaHypermutationIds)]
gliomaNONHypermutatorMaf = mafAdj[mafAdj['Tumor_Sample_Barcode'].isin(gliomaNotHypermutationIds)]
occurenceCounterHypermutators, rankingDictHypermutators, fractionalDictHypermutators = maf_analysis_utils.enumerate_top_n_oncogenic_mutated_genes_across_cohort(gliomaHypermutatorMaf, n=None)
occurenceCounterNONHypermutators, rankingDictNONHypermutators, fractionalDictNONHypermutators = maf_analysis_utils.enumerate_top_n_oncogenic_mutated_genes_across_cohort(gliomaNONHypermutatorMaf, n=None)

gliomaGeneDf = make_gene_based_df(fractionalDictHypermutators, occurenceCounterHypermutators, fractionalDictNONHypermutators, occurenceCounterNONHypermutators, rankingDictHypermutators, cdsGeneLengthDict)
gliomaGeneDf.to_csv('~/Desktop/WORK/dataForLocalPlotting/gliomaHypermutationComparisson.tsv', sep='\t', index=False)

#MELANOMA
melanomaSigs = impactSigs[impactSigs['cancer_type'] == 'Melanoma']
melanomaHypermutationIds = set(melanomaSigs[melanomaSigs['Nmut_Mb'] > hypermutationNmutMbThresh]['Tumor_Sample_Barcode'])
melanomaNotHypermutationIds = set(melanomaSigs[melanomaSigs['Nmut_Mb'] <= hypermutationNmutMbThresh]['Tumor_Sample_Barcode'])
melanomaHypermutatorMaf = mafAdj[mafAdj['Tumor_Sample_Barcode'].isin(melanomaHypermutationIds)]
melanomaNONHypermutatorMaf = mafAdj[mafAdj['Tumor_Sample_Barcode'].isin(melanomaNotHypermutationIds)]
occurenceCounterHypermutators, rankingDictHypermutators, fractionalDictHypermutators = maf_analysis_utils.enumerate_top_n_oncogenic_mutated_genes_across_cohort(melanomaHypermutatorMaf, n=None)
occurenceCounterNONHypermutators, rankingDictNONHypermutators, fractionalDictNONHypermutators = maf_analysis_utils.enumerate_top_n_oncogenic_mutated_genes_across_cohort(melanomaNONHypermutatorMaf, n=None)

melanomaGeneDf = make_gene_based_df(fractionalDictHypermutators, occurenceCounterHypermutators, fractionalDictNONHypermutators, occurenceCounterNONHypermutators, rankingDictHypermutators, cdsGeneLengthDict)
melanomaGeneDf.to_csv('~/Desktop/WORK/dataForLocalPlotting/melanomaHypermutationComparisson.tsv', sep='\t', index=False)

#BLADDER
bladderSigs = impactSigs[impactSigs['cancer_type'] == 'Bladder Cancer']
bladderHypermutationIds = set(bladderSigs[bladderSigs['Nmut_Mb'] > hypermutationNmutMbThresh]['Tumor_Sample_Barcode'])
bladderNotHypermutationIds = set(bladderSigs[bladderSigs['Nmut_Mb'] <= hypermutationNmutMbThresh]['Tumor_Sample_Barcode'])
bladderHypermutatorMaf = mafAdj[mafAdj['Tumor_Sample_Barcode'].isin(bladderHypermutationIds)]
bladderNONHypermutatorMaf = mafAdj[mafAdj['Tumor_Sample_Barcode'].isin(bladderNotHypermutationIds)]
occurenceCounterHypermutators, rankingDictHypermutators, fractionalDictHypermutators = maf_analysis_utils.enumerate_top_n_oncogenic_mutated_genes_across_cohort(bladderHypermutatorMaf, n=None)
occurenceCounterNONHypermutators, rankingDictNONHypermutators, fractionalDictNONHypermutators = maf_analysis_utils.enumerate_top_n_oncogenic_mutated_genes_across_cohort(bladderNONHypermutatorMaf, n=None)

bladderGeneDf = make_gene_based_df(fractionalDictHypermutators, occurenceCounterHypermutators, fractionalDictNONHypermutators, occurenceCounterNONHypermutators, rankingDictHypermutators, cdsGeneLengthDict)
bladderGeneDf.to_csv('~/Desktop/WORK/dataForLocalPlotting/bladderHypermutationComparisson.tsv', sep='\t', index=False)



len(bladderHypermutationIds), len(bladderNotHypermutationIds)


