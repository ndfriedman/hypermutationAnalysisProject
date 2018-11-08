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
#We also need to count the number of mutations that occur at NOTaging and NOTtmz motifs
allTrinucs = analysis_utils.get_all_possible_quadNucs()
nonTmzNonAgingTrinucs = allTrinucs - spectraEnrichmentDict['Signature.11'] - spectraEnrichmentDict['Signature.1']
gliomaSummaryInfoOther = maf_analysis_utils.summarize_signature_attribution_for_case(gliomaAnalysisMuts, nonTmzNonAgingTrinucs)

gliomaSummaryInfoTMZ['mutSource'] = 'TMZ'
gliomaSummaryInfoAging['mutSource'] = 'Aging'
gliomaSummaryInfoOther['mutSource'] = 'Other'

concatDf = pd.concat([gliomaSummaryInfoTMZ, gliomaSummaryInfoAging, gliomaSummaryInfoOther])

#LETS plot this shit!
concatDf.to_csv('')


