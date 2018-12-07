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

#Focus on bladder for now

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
endometrialSigs = impactSigs[impactSigs['cancer_type'] == 'Endometrial Cancer'] #Ignore MMR & POLE bladders

endometrialIds = set(endometrialSigs['Tumor_Sample_Barcode']) 

endometrialAnalysisMuts = mafWithInfo[mafWithInfo['Tumor_Sample_Barcode'].isin(endometrialIds)]
endometrialAnalysisMuts['quadNuc'] = endometrialAnalysisMuts.apply(lambda row: mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)

sig10Spectra = spectraEnrichmentDict['Signature.10']
sig14Spectra = spectraEnrichmentDict['Signature.14']
mmrSpectra = spectraEnrichmentDict['Signature.MMR']
agingSpectra = spectraEnrichmentDict['Signature.1']

allTrinucs = analysis_utils.get_all_possible_quadNucs()
agingAndPoleSpectra = sig10Spectra & agingSpectra
poleOnlySpectra = sig10Spectra - agingAndPoleSpectra
sig14Spectra = sig14Spectra - poleOnlySpectra #Remove the one enriched POLE signature from the sig14 spectra
mmrORSig14Spectra = sig14Spectra | mmrSpectra
mmrAndAgingSpectra = mmrORSig14Spectra & agingSpectra
mmrOnlySpectra = mmrORSig14Spectra - agingSpectra
agingOnlySpectra = agingSpectra - mmrAndAgingSpectra - agingAndPoleSpectra
allOtherSpectra = allTrinucs - poleOnlySpectra - agingOnlySpectra - mmrOnlySpectra - mmrAndAgingSpectra - agingAndPoleSpectra

len(allOtherSpectra), len(agingOnlySpectra), len(poleOnlySpectra), len(mmrAndAgingSpectra), len(mmrOnlySpectra), len(agingAndPoleSpectra)

#There are five categories: POLE only spectra (2), MMR/14 only spectra (12), MMR&Aging Spectra (3), pole&agingSpectra(1) and All other spectra (78) 

endometrialSummaryInfoPOLE = maf_analysis_utils.summarize_signature_attribution_for_case(endometrialAnalysisMuts, poleOnlySpectra)
endometrialSummaryInfoMMR = maf_analysis_utils.summarize_signature_attribution_for_case(endometrialAnalysisMuts, mmrOnlySpectra)
endometrialSummaryInfoMMRAging = maf_analysis_utils.summarize_signature_attribution_for_case(endometrialAnalysisMuts, mmrAndAgingSpectra)
endometrialSummaryInfoPOLEAging = maf_analysis_utils.summarize_signature_attribution_for_case(endometrialAnalysisMuts, agingAndPoleSpectra)
endometrialSummaryInfoOther = maf_analysis_utils.summarize_signature_attribution_for_case(endometrialAnalysisMuts, allOtherSpectra)

#Rename columns for plotting
endometrialSummaryInfoPOLE['mutSource'] = 'POLE'
endometrialSummaryInfoMMR['mutSource'] = 'MMR'
endometrialSummaryInfoMMRAging['mutSource'] = 'MMR or Aging'
endometrialSummaryInfoPOLEAging['mutSource'] = 'POLE or Aging'
endometrialSummaryInfoOther['mutSource'] = 'Other'

concatDf = pd.concat([endometrialSummaryInfoPOLE, endometrialSummaryInfoMMR, endometrialSummaryInfoMMRAging, endometrialSummaryInfoPOLEAging, endometrialSummaryInfoOther])

#LETS plot this shit!
concatDf.to_csv('~/Desktop/dataForLocalPlotting/endometrialMutAttributionData.tsv', sep='\t', index=False)














