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

#maf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/annotatedOncoPlusHotspotMafAllImpact_trinuc')
#maf = maf_analysis_utils.add_mut_effect_summary_col(maf)
#impactSigs = pd.read_table(pathPrefix + '/ifs/res/taylorlab/impact_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
#impactSigsInMaf = impactSigs[impactSigs['Tumor_Sample_Barcode'].isin(set(maf['Tumor_Sample_Barcode']))]
#overlapIds = set((impactSigsInMaf['Tumor_Sample_Barcode']))

#sigsOverlap = impactSigs[impactSigs['Tumor_Sample_Barcode'].isin(overlapIds)]
#mafOverlap = maf[maf['Tumor_Sample_Barcode'].isin(overlapIds)]
#sigsOverlap.to_csv(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/hypermutationAnalysisProjectData/hypermutationAnalysisProjectSignatures.tsv', sep='\t', index=False)
#mafOverlap.to_csv(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/hypermutationAnalysisProjectData/hypermutationAnalysisProjectMaf.tsv', sep='\t', index=False)

maf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/hypermutationAnalysisProjectData/hypermutationAnalysisProjectMaf.tsv')
sigs = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/hypermutationAnalysisProjectData/hypermutationAnalysisProjectSignatures.tsv')

oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])

np.nanmean(sigs['Nmut_Mb'])

hypermutationThreshold = 80 #the threshold for cases to be considered hypermutators
hypermutatedIds = set(sigs[sigs['Nmut_Mb'] > hypermutationThreshold]['Tumor_Sample_Barcode'])
hypermutatedMaf = maf[maf['Tumor_Sample_Barcode'].isin(hypermutatedIds)]
nonHypermutatedMaf = maf[~maf['Tumor_Sample_Barcode'].isin(hypermutatedIds)]

maf.shape
nonHypermutatedMaf.shape
len(set(nonHypermutatedMaf['Tumor_Sample_Barcode']))

mafAnalyze = nonHypermutatedMaf

nMafMuts = mafAnalyze.shape[0]
nMafCases = len(set(mafAnalyze['Tumor_Sample_Barcode']))
hotspotMutMaf = mafAnalyze[mafAnalyze['is-a-hotspot'] == 'Y']
print 'hotspot', hotspotMutMaf.shape[0], 1.0*hotspotMutMaf.shape[0]/nMafMuts, len(set(hotspotMutMaf['Tumor_Sample_Barcode'])), 1.0*len(set(hotspotMutMaf['Tumor_Sample_Barcode']))/nMafCases

oncogenicMutMaf = mafAnalyze[mafAnalyze['oncogenic'].isin(oncogenicMutColNames)]
print 'oncogenic', oncogenicMutMaf.shape[0], 1.0*oncogenicMutMaf.shape[0]/nMafMuts, len(set(oncogenicMutMaf['Tumor_Sample_Barcode'])), 1.0*len(set(oncogenicMutMaf['Tumor_Sample_Barcode']))/nMafCases

consequenceMutMaf = mafAnalyze[(mafAnalyze['Consequence'] == 'stop_gained') | (mafAnalyze['Consequence'] == 'frameshift_variant')]
print 'consequential', consequenceMutMaf.shape[0], 1.0*consequenceMutMaf.shape[0]/nMafMuts, len(set(consequenceMutMaf['Tumor_Sample_Barcode'])), 1.0*len(set(consequenceMutMaf['Tumor_Sample_Barcode']))/nMafCases

knownOncTruncMaf = mafAnalyze[mafAnalyze['mutKnowOncogenicOrTruncating'] == True]
print 'knownOncTrunc', knownOncTruncMaf.shape[0], 1.0*knownOncTruncMaf.shape[0]/nMafMuts, len(set(knownOncTruncMaf['Tumor_Sample_Barcode'])), 1.0*len(set(knownOncTruncMaf['Tumor_Sample_Barcode']))/nMafCases

print len(set(hypermutatedMaf['HGVSp']))
print len(set(hypermutatedMaf[hypermutatedMaf['is-a-hotspot'] == 'Y']['HGVSp']))
print len(set(hypermutatedMaf[hypermutatedMaf['oncogenic'].isin(oncogenicMutColNames)]['HGVSp']))
print len(set(hypermutatedMaf[(hypermutatedMaf['Consequence'] == 'stop_gained') | (hypermutatedMaf['Consequence'] == 'frameshift_variant')]['HGVSp']))
print len(set(hypermutatedMaf[hypermutatedMaf['mutKnowOncogenicOrTruncating'] == True]['HGVSp']))


