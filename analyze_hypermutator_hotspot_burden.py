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

#returns a dataframe mapping case names, n oncogenic mutations and n hotspots, fraction hotspots at signature enriched motif per case
def enumerate_case_mutation_info_summary(df, enrichedSigMotifs):
    listOfDicts = []
    cases = set(df['Tumor_Sample_Barcode'])
    for case in cases:
        localD = {}
        caseDf = df[df['Tumor_Sample_Barcode'] == case]
        
        #TODO EXTEND THIS ANALYSIS TO INCORP FRAC HOTSPOTS at enriched motif and FRAC oncogenic muts at enriched MOTIF
        
        oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic']) #enumerate col names for likely oncogenic mutations
        
        hotspotMutDf = caseDf[caseDf['is-a-hotspot'] == 'Y']
        oncogenicMutDf = caseDf[caseDf['oncogenic'].isin(oncogenicMutColNames)]
        oncogenicOrHotspotMutations = caseDf[(caseDf['oncogenic'].isin(oncogenicMutColNames)) | (caseDf['is-a-hotspot'] == 'Y')]
        
        nHotspotMutations = hotspotMutDf.shape[0]
        nOncogenicMutations = oncogenicMutDf.shape[0]
        nOncogenicOrHotspotMutations = oncogenicOrHotspotMutations.shape[0] #we need to count this separately because they may overlap 
        
        fracOncogenicMutationsAtEnrichedMotif = None
        fracHotpsotMutationsAtEnrichedMotif = None
        fracDriverMutationsAtEnrichedMotif = None
        if nHotspotMutations > 0:
            fracHotpsotMutationsAtEnrichedMotif = 1.0*hotspotMutDf[hotspotMutDf['quadNuc'].isin(enrichedSigMotifs)].shape[0]/nHotspotMutations
        if nOncogenicOrHotspotMutations > 0:
            fracOncogenicMutationsAtEnrichedMotif = 1.0*oncogenicMutDf[oncogenicMutDf['quadNuc'].isin(enrichedSigMotifs)].shape[0]/nOncogenicMutations
        if nOncogenicOrHotspotMutations > 0:
            fracDriverMutationsAtEnrichedMotif = 1.0*oncogenicOrHotspotMutations[oncogenicOrHotspotMutations['quadNuc'].isin(enrichedSigMotifs)].shape[0]/nOncogenicOrHotspotMutations
        
        nMut = caseDf.shape[0]
        nMutToHotspotRatio = None
        nMutToOncogenicRatio = None
        nMutToOncogenicAndHotspotRatio = None
        if nHotspotMutations > 0:
            nMutToHotspotRatio = 1.0*nMut/nHotspotMutations
        if nOncogenicMutations > 0:
            nMutToOncogenicRatio = 1.0*nMut/nOncogenicMutations
        if nOncogenicOrHotspotMutations > 0:
            nMutToOncogenicAndHotspotRatio = 1.0*nMut/nOncogenicOrHotspotMutations
        
        
        #add in all the information to the local dict
        localD['Tumor_Sample_Barcode'] = case
        localD['nHotspots'] = nHotspotMutations
        localD['nOncogenicMutations'] = nOncogenicMutations
        localD['nOncogenicOrHotspotMutations'] = nOncogenicOrHotspotMutations
        localD['fracOncogenicMutationsAtEnrichedMotif'] = fracOncogenicMutationsAtEnrichedMotif
        localD['fracHotspotMutationsAtEnrichedMotif'] = fracHotpsotMutationsAtEnrichedMotif
        localD['fracDriverMutationsAtEnrichedMotif'] = fracDriverMutationsAtEnrichedMotif
        localD['nMutToHotspotRatio'] = nMutToHotspotRatio
        localD['nMutToOncogenicRatio'] = nMutToOncogenicRatio
        localD['nMutToOncogenicAndHotspotRatio'] = nMutToOncogenicAndHotspotRatio
        localD['Nmut'] = nMut
        
        listOfDicts.append(localD)
      
    df = pd.DataFrame(listOfDicts)
    return df

mafWithInfo = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/annotatedOncoPlusHotspotMafAllImpact_trinuc')
impactSigs = pd.read_table(pathPrefix + '/ifs/res/taylorlab/impact_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
impactSigs['pid'] = impactSigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
impactSigs['cancer_type'] = impactSigs['pid'].apply(lambda x: cDict[x] if x in cDict else None)

#adjust column names to make the 
renameDict = {key:value for (key,value) in [('mean_' + str(i), 'Signature.' + str(i)) for i in range(1,31)]}
impactSigs = impactSigs.rename(columns=renameDict)
impactSigs = mutationSigUtils.merge_signature_columns(impactSigs, mode='Stratton', drop=False, smokingMerge=True, confidence=False, mean=True, prefix='Signature.')

#get a dictionary that maps the trinuceleotide contexts at which certain signatures are enriched
spectraEnrichmentDict = mutationSigUtils.get_enriched_spectra_for_signatures(spectraSignificanceThresh=.05, pathPrefix='/Users/friedman/Desktop/mnt',
	signaturesToIgnore= #ignore signatures we dont care about 
	set(['Signature.5','Signature.8','Signature.9','Signature.12','Signature.16','Signature.19','Signature.22','Signature.23','Signature.24','Signature.25','Signature.27','Signature.28','Signature.29','Signature.30']))

#add quad nuc info
mafWithInfo['quadNuc'] = mafWithInfo.apply(lambda row: mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)

hypermutationThreshold = 75
signatureDetectionThreshold = .25
#make sets of the ids of hyper mutators 

hypermutatorPoleCases = set(impactSigs[(impactSigs['Signature.10'] >= signatureDetectionThreshold) & (impactSigs['Nmut'] >= hypermutationThreshold)]['Tumor_Sample_Barcode'])
hypermutatorMMRCases = set(impactSigs[(impactSigs['Signature.MMR'] >= signatureDetectionThreshold) & (impactSigs['Nmut'] >= hypermutationThreshold)]['Tumor_Sample_Barcode'])
hypermutatorUVCases = set(impactSigs[(impactSigs['Signature.7'] >= signatureDetectionThreshold) & (impactSigs['Nmut'] >= hypermutationThreshold)]['Tumor_Sample_Barcode'])
hypermutatorAPOBECCases = set(impactSigs[(impactSigs['Signature.APOBEC'] >= signatureDetectionThreshold) & (impactSigs['Nmut'] >= hypermutationThreshold)]['Tumor_Sample_Barcode'])
hypermutatorTMZCases = set(impactSigs[(impactSigs['Signature.11'] >= signatureDetectionThreshold) & (impactSigs['Nmut'] >= hypermutationThreshold)]['Tumor_Sample_Barcode'])
hypermutatorMMRCases = hypermutatorMMRCases - hypermutatorPoleCases #(cases with MMR plus POLE are classified as just MMR)

#make sets of controls
#endometrial MSS: endometrials with < 75 muts and <25% mmr and msi signature
endometrialMSSNonPole = set(impactSigs[(impactSigs['cancer_type'] == 'Endometrial Cancer') & ((impactSigs['Signature.10'] < signatureDetectionThreshold) & (impactSigs['Nmut'] < hypermutationThreshold) & (impactSigs['Signature.MMR'] < signatureDetectionThreshold))]['Tumor_Sample_Barcode'])
#colorectal MMS: colorectal cancers with <75 muts and <25% mmr and msi signature
colorectalMSS = set(impactSigs[(impactSigs['cancer_type'] == 'Colorectal Cancer')&(impactSigs['Signature.MMR'] < signatureDetectionThreshold)&(impactSigs['Nmut'] < hypermutationThreshold)]['Tumor_Sample_Barcode'])
#non UV melanoma
nonUVMelanoma = set(impactSigs[(impactSigs['cancer_type'] == 'Melanoma')&(impactSigs['Signature.7'] < signatureDetectionThreshold)&(impactSigs['Nmut'] < hypermutationThreshold)]['Tumor_Sample_Barcode'])
#nonApobecBladder (bc bladder is the most common apobec hypermutator)
nonApobecBladder = set(impactSigs[(impactSigs['cancer_type'] == 'Bladder Cancer')&(impactSigs['Signature.APOBEC'] < signatureDetectionThreshold)&(impactSigs['Nmut'] < hypermutationThreshold)]['Tumor_Sample_Barcode'])
#non TMZ glioma
nonTMZGlioma = set(impactSigs[(impactSigs['cancer_type'] == 'Glioma')&(impactSigs['Signature.11'] < signatureDetectionThreshold)&(impactSigs['Nmut'] < hypermutationThreshold)]['Tumor_Sample_Barcode'])

#Optionally make a summary DF for everything
mafWithInfo.shape
fullSummaryDf.shape = enumerate_case_mutation_info_summary(mafWithInfo, {})

#Make summary DFs for all the cases we care about
poleSummaryDf = enumerate_case_mutation_info_summary(mafWithInfo[mafWithInfo['Tumor_Sample_Barcode'].isin(hypermutatorPoleCases)], spectraEnrichmentDict['Signature.10'])
mmrSummaryDf = enumerate_case_mutation_info_summary(mafWithInfo[mafWithInfo['Tumor_Sample_Barcode'].isin(hypermutatorMMRCases)], spectraEnrichmentDict['Signature.MMR'])
uvSummaryDf = enumerate_case_mutation_info_summary(mafWithInfo[mafWithInfo['Tumor_Sample_Barcode'].isin(hypermutatorUVCases)], spectraEnrichmentDict['Signature.7'])
apobecSummaryDf = enumerate_case_mutation_info_summary(mafWithInfo[mafWithInfo['Tumor_Sample_Barcode'].isin(hypermutatorAPOBECCases)], spectraEnrichmentDict['Signature.APOBEC'])
tmzSummaryDf = enumerate_case_mutation_info_summary(mafWithInfo[mafWithInfo['Tumor_Sample_Barcode'].isin(hypermutatorTMZCases)], spectraEnrichmentDict['Signature.11'])
endometrialMSSNonPoleSummaryDf = enumerate_case_mutation_info_summary(mafWithInfo[mafWithInfo['Tumor_Sample_Barcode'].isin(endometrialMSSNonPole)], spectraEnrichmentDict['Signature.1'])
colorectalMSSSummaryDf = enumerate_case_mutation_info_summary(mafWithInfo[mafWithInfo['Tumor_Sample_Barcode'].isin(colorectalMSS)], spectraEnrichmentDict['Signature.1'])
nonUVMelanomaSummaryDf = enumerate_case_mutation_info_summary(mafWithInfo[mafWithInfo['Tumor_Sample_Barcode'].isin(nonUVMelanoma)], spectraEnrichmentDict['Signature.1'])
nonApobecBladderSummaryDf = enumerate_case_mutation_info_summary(mafWithInfo[mafWithInfo['Tumor_Sample_Barcode'].isin(nonApobecBladder)], spectraEnrichmentDict['Signature.1'])
nonTMZGliomaSummaryDf = enumerate_case_mutation_info_summary(mafWithInfo[mafWithInfo['Tumor_Sample_Barcode'].isin(nonTMZGlioma)], spectraEnrichmentDict['Signature.1'])

#SET ALL the labels for dataframes before concatenating them Also set numbers for plotOrdering
#higher numbers go first in R
poleSummaryDf['label'] = 'POLE'
poleSummaryDf['plotOrdering'] = 10
mmrSummaryDf['label'] = 'MMR'
mmrSummaryDf['plotOrdering'] = 8
uvSummaryDf['label'] = 'UV'
uvSummaryDf['plotOrdering'] = 4
apobecSummaryDf['label'] = 'APOBEC'
apobecSummaryDf['plotOrdering'] = 2
tmzSummaryDf['label'] = 'TMZ'
tmzSummaryDf['plotOrdering'] = 6
endometrialMSSNonPoleSummaryDf['label'] = 'EndometrialMSS'
endometrialMSSNonPoleSummaryDf['plotOrdering'] = 9
colorectalMSSSummaryDf['label'] = 'ColorectalMSS'
colorectalMSSSummaryDf['plotOrdering'] = 7
nonUVMelanomaSummaryDf['label'] = 'nonUVMelanoma'
nonUVMelanomaSummaryDf['plotOrdering'] = 3
nonApobecBladderSummaryDf['label'] = 'nonApobecBladder'
nonApobecBladderSummaryDf['plotOrdering'] = 1
nonTMZGliomaSummaryDf['label'] = 'nonTMZGlioma'
nonTMZGliomaSummaryDf['plotOrdering'] = 5

concatedDf = pd.concat([poleSummaryDf, mmrSummaryDf, uvSummaryDf, apobecSummaryDf, tmzSummaryDf, endometrialMSSNonPoleSummaryDf, colorectalMSSSummaryDf, nonUVMelanomaSummaryDf, nonApobecBladderSummaryDf, nonTMZGliomaSummaryDf])
concatedDf.to_csv('~/Desktop/dataForLocalPlotting/mutburdenBoxplot.tsv', sep='\t', index=False)




hypermutatorIds = set(impactSigs[impactSigs['Nmut'] >= 75]['pid'])
len(hypermutatorIds)

for v in list(tmzSummaryDf['Tumor_Sample_Barcode']):
    print v

#TODO make the endometrials MSS 
#add MSI, UV, and TMZ
#add frac hotspots attributable to dominant signature????

Counter(mafWithInfo[(mafWithInfo['Tumor_Sample_Barcode'].isin(hypermutatorPoleCases)) &(mafWithInfo['is-a-hotspot'] == 'Y')]['Hugo_Symbol']).most_common(10)

Counter(mafWithInfo[(mafWithInfo['Tumor_Sample_Barcode'].isin(endometrialMSSNonPole)) &(mafWithInfo['is-a-hotspot'] == 'Y') & (mafWithInfo['Hugo_Symbol'] == 'PIK3CA')]['amino_acid_change'])

mafWithInfo.columns.values

for case in endometrialMSSNonPole:
    print case

mafWithInfoPoleHotspots = mafWithInfoPole[mafWithInfoPole['is-a-hotspot'] == 'Y']
enumerate_n_hotspots_per_case(mafWithInfoPoleHotspots)

