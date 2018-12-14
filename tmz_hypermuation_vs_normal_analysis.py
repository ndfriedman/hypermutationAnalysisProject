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

mafWithInfo = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/annotatedOncoPlusHotspotMafAllImpact_trinuc')
impactSigs = pd.read_table(pathPrefix + '/ifs/res/taylorlab/impact_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
impactSigs['pid'] = impactSigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
impactSigs['cancer_type'] = impactSigs['pid'].apply(lambda x: cDict[x] if x in cDict else None)

#adjust column names to make the 
renameDict = {key:value for (key,value) in [('mean_' + str(i), 'Signature.' + str(i)) for i in range(1,31)]}
impactSigs = impactSigs.rename(columns=renameDict)
impactSigs = mutationSigUtils.merge_signature_columns(impactSigs, mode='Stratton', drop=False, smokingMerge=True, confidence=False, mean=True, prefix='Signature.')

spectraEnrichmentDict = mutationSigUtils.get_enriched_spectra_for_signatures(spectraSignificanceThresh=.05, pathPrefix='/Users/friedman/Desktop/mnt',
	signaturesToIgnore= #ignore signatures we dont care about 
	set(['Signature.5','Signature.8','Signature.9','Signature.12','Signature.16','Signature.19','Signature.22','Signature.23','Signature.24','Signature.25','Signature.27','Signature.28','Signature.29','Signature.30']))

gliomaSigs = impactSigs[impactSigs['cancer_type'] == 'Glioma']
gliomaIds = set(gliomaSigs['Tumor_Sample_Barcode'])
#be a bit finnicky with the analysis: if there are more than 100 muts present they have to be TMZ hypermutated
gliomaAnalysisIds = set(gliomaSigs[(gliomaSigs['Nmut'] < 100) | ((gliomaSigs['Nmut'] >= 100) & (gliomaSigs['Signature.11'] > .25))]['Tumor_Sample_Barcode'])
gliomaAnalysisMuts = mafWithInfo[mafWithInfo['Tumor_Sample_Barcode'].isin(gliomaAnalysisIds)]
gliomaAnalysisMuts['quadNuc'] = gliomaAnalysisMuts.apply(lambda row: mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)

gliomaSummaryInfoTMZ = maf_analysis_utils.summarize_signature_attribution_for_case(gliomaAnalysisMuts, spectraEnrichmentDict['Signature.11'])
gliomaSummaryInfoAging = maf_analysis_utils.summarize_signature_attribution_for_case(gliomaAnalysisMuts, spectraEnrichmentDict['Signature.1'])
#we calculate the portion of TMZ mutations at the weakly enriched TMZ motif
weakTMZmotif = set(['ACTA', 'CCTA', 'GCTA', 'TCTA'])
gliomaSummaryInfoWeakTMZ = maf_analysis_utils.summarize_signature_attribution_for_case(gliomaAnalysisMuts, weakTMZmotif)
#We also need to count the number of mutations that occur at NOTaging and NOTtmz motifs
allTrinucs = analysis_utils.get_all_possible_quadNucs()
nonTmzNonAgingTrinucs = allTrinucs - spectraEnrichmentDict['Signature.11'] - spectraEnrichmentDict['Signature.1'] - weakTMZmotif
gliomaSummaryInfoOther = maf_analysis_utils.summarize_signature_attribution_for_case(gliomaAnalysisMuts, nonTmzNonAgingTrinucs)

gliomaSummaryInfoTMZ['mutSource'] = 'TMZ (top 8 tmz most favored trinucs)'
gliomaSummaryInfoWeakTMZ['mutSource'] = 'probably TMZ (TMZ 9th-12th most favored trinucs)'
gliomaSummaryInfoAging['mutSource'] = 'Aging (aging top four favored trinucs)'
gliomaSummaryInfoOther['mutSource'] = 'Other (non C>T mutations)'

concatDf = pd.concat([gliomaSummaryInfoTMZ, gliomaSummaryInfoAging, gliomaSummaryInfoOther, gliomaSummaryInfoWeakTMZ])

#LETS plot this shit!
concatDf.to_csv('~/Desktop/dataForLocalPlotting/gliomaMutAttributionData.tsv', sep='\t', index=False)




###############################################################################
#Do gene comparissons
tmzStrongMotif = spectraEnrichmentDict['Signature.11']
weakTMZmotif = set(['ACTA', 'CCTA', 'GCTA', 'TCTA'])

gliomaNonTMZhyperIDs = set(gliomaSigs[(gliomaSigs['Signature.11'] < .25) & (gliomaSigs['Nmut'] < 100)]['Tumor_Sample_Barcode'])
gliomaTMZhyperIDs = set(gliomaSigs[(gliomaSigs['Signature.11'] >= .25) & (gliomaSigs['Nmut'] >= 100)]['Tumor_Sample_Barcode'])

for tid in gliomaTMZhyperIDs:
    print tid

oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
oncogenicGliomaMuts = gliomaAnalysisMuts[gliomaAnalysisMuts['oncogenic'].isin(oncogenicMutColNames)]
nonOncogenicGliomaMuts = gliomaAnalysisMuts[~gliomaAnalysisMuts['oncogenic'].isin(oncogenicMutColNames)]

nTMZHypermutated = len(gliomaTMZhyperIDs)
oncogenicMutsInTMZHypers = oncogenicGliomaMuts[oncogenicGliomaMuts['Tumor_Sample_Barcode'].isin(gliomaTMZhyperIDs)]
oncogenicMutsInNonHyperGlioma = oncogenicGliomaMuts[oncogenicGliomaMuts['Tumor_Sample_Barcode'].isin(gliomaNonTMZhyperIDs)]
top50MostMutatedGenesInHyper = [i[0] for i in Counter(oncogenicMutsInTMZHypers['Hugo_Symbol']).most_common(50)]

#calculate fractions mutated for the non hypermuated cases
nNonHyperGlioma = len(gliomaNonTMZhyperIDs)
nonHyperFracDict = dict()
for gene in top50MostMutatedGenesInHyper:
    nCasesWithGeneMutated = len(set(oncogenicMutsInNonHyperGlioma[oncogenicMutsInNonHyperGlioma['Hugo_Symbol'] == gene]['Tumor_Sample_Barcode']))
    nonHyperFracDict[gene] = 1.0*nCasesWithGeneMutated/nNonHyperGlioma


geneLengthDict = analysis_utils.get_gene_length_info(bedFilePath = pathPrefix + '/ifs/res/pwg/data/gencode/gencode.v19.all_gene_bounds.bed')

concatDfList = []
for gene in top50MostMutatedGenesInHyper:
    curSumDf = maf_analysis_utils.asses_per_case_mut_info_for_gene(oncogenicGliomaMuts[oncogenicGliomaMuts['Tumor_Sample_Barcode'].isin(gliomaTMZhyperIDs)], gene, tmzStrongMotif | weakTMZmotif)
    curSumDf['ordering'] = curSumDf.shape[0]
    if gene in geneLengthDict:
        curSumDf['geneLength'] = geneLengthDict[gene]
    else:
        curSumDf['geneLength'] = None
    fracMutatedInHypermutators = 1.0*curSumDf.shape[0]/nTMZHypermutated
    #if nonHyperFracDict[gene] == 0:
    #curSumDf['ratio'] = None
    #else:
    curSumDf['ratio'] = nonHyperFracDict[gene]/fracMutatedInHypermutators
    concatDfList.append(curSumDf)
concatDf = pd.concat(concatDfList)

print concatDf['ratio']

concatDf.to_csv('~/Desktop/dataForLocalPlotting/gliomaHypermutationDistribution.tsv', sep='\t', index=False)

###########################################

print analysis_utils.normalize_counter(Counter(gliomaAnalysisMuts[gliomaAnalysisMuts['Tumor_Sample_Barcode'].isin(gliomaNonTMZhyperIDs)]['Hugo_Symbol']), nDigitsRound=3).most_common(20)
print '_________'
print analysis_utils.normalize_counter(Counter(gliomaAnalysisMuts[gliomaAnalysisMuts['Tumor_Sample_Barcode'].isin(gliomaTMZhyperIDs)]['Hugo_Symbol']), nDigitsRound=3).most_common(20)


print analysis_utils.normalize_counter(Counter(oncogenicGliomaMuts[oncogenicGliomaMuts['Tumor_Sample_Barcode'].isin(gliomaNonTMZhyperIDs)]['Hugo_Symbol']), nDigitsRound=3).most_common(20)
print '_________'
print analysis_utils.normalize_counter(Counter(oncogenicGliomaMuts[oncogenicGliomaMuts['Tumor_Sample_Barcode'].isin(gliomaTMZhyperIDs)]['Hugo_Symbol']), nDigitsRound=3).most_common(20)

#TODO properly get per case frequencies for each gene with dropping duplicates
#print analysis_utils.normalize_counter(Counter(oncogenicGliomaMuts[oncogenicGliomaMuts['Tumor_Sample_Barcode'].isin(gliomaNonTMZhyperIDs)]['Hugo_Symbol'].drop_duplicates(subset=['Tumor_Sample_Barcode'])['Hugo_Symbol']), nDigitsRound=3).most_common(20)
#print '_________'
#print analysis_utils.normalize_counter(Counter(oncogenicGliomaMuts[oncogenicGliomaMuts['Tumor_Sample_Barcode'].isin(gliomaTMZhyperIDs)].drop_duplicates(subset=['Tumor_Sample_Barcode'])['Hugo_Symbol']), nDigitsRound=3).most_common(20)

oncogenicGliomaMuts[oncogenicGliomaMuts['Tumor_Sample_Barcode'].isin(gliomaNonTMZhyperIDs)].drop_duplicates(subset=['Tumor_Sample_Barcode']).shape
v = oncogenicGliomaMuts[oncogenicGliomaMuts['Tumor_Sample_Barcode'].isin(gliomaTMZhyperIDs)]
print v[v['Hugo_Symbol'] == 'TP53'].shape
print len(set(v[v['Hugo_Symbol'] == 'TP53']['Tumor_Sample_Barcode']))
d = v[v['Hugo_Symbol'] == 'TP53'].drop_duplicates(subset=['Tumor_Sample_Barcode'])
print d[d['Hugo_Symbol'] == 'TP53'].shape

for x in gliomaTMZhyperIDs:
    print x

ptenNonOncogenicTMZhyper = ptenNonOncogenicGliomaMuts[ptenNonOncogenicGliomaMuts['Tumor_Sample_Barcode'].isin(gliomaTMZhyperIDs)]
ptenOncogenicTMZhyper = ptenOncogenicGliomaMuts[ptenOncogenicGliomaMuts['Tumor_Sample_Barcode'].isin(gliomaTMZhyperIDs)]

ptenNonOncogenicNonTMZ = ptenNonOncogenicGliomaMuts[ptenNonOncogenicGliomaMuts['Tumor_Sample_Barcode'].isin(gliomaNonTMZhyperIDs)]


print ptenOncogenicTMZhyper[ptenOncogenicTMZhyper['quadNuc'].isin(tmzStrongMotif)].shape, ptenOncogenicTMZhyper.shape

#hotspotMaf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/hotspotReducedAnalysis10-19.tsv')
hotspotMaf['HGVSp']




































hotspotMaf.columns.values













