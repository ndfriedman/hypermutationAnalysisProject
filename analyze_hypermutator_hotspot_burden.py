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
    cntr = 0
    for case in cases:
        if cntr%100 == 0: 
            print cntr, len(cases)
        cntr += 1
        
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
mafWithInfo['pid'] = mafWithInfo['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
mafWithInfo['cancer_type'] = mafWithInfo['pid'].apply(lambda x: cDict[x] if x in cDict else None)

#adjust column names to make the 
renameDict = {key:value for (key,value) in [('mean_' + str(i), 'Signature.' + str(i)) for i in range(1,31)]}
impactSigs = impactSigs.rename(columns=renameDict)
impactSigs = mutationSigUtils.merge_signature_columns(impactSigs, mode='Stratton', drop=False, smokingMerge=True, confidence=False, mean=True, prefix='Signature.')

nmut_MbDict = dict(zip(impactSigs['Tumor_Sample_Barcode'], impactSigs['Nmut_Mb']))

cancerTypesToFocusOn = set(['Colorectal Cancer', 'Non-Small Cell Lung Cancer', 'Glioma', 'Melanoma', 'Endometrial Cancer', 'Bladder Cancer'])

#mafLimitedToCancerType = mafWithInfo[mafWithInfo['cancer_type'].isin(cancerTypesToFocusOn)]

#ADD IN INFORMATION ABOUT QUADNUCs
#mafLimitedToCancerType['quadNuc'] = mafLimitedToCancerType.apply(lambda row: mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)
mafWithInfo['quadNuc'] = mafWithInfo.apply(lambda row: mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)

#infoDf  = enumerate_case_mutation_info_summary(mafLimitedToCancerType, set([]))
infoDf  = enumerate_case_mutation_info_summary(mafWithInfo, set([]))

infoDf['Nmut_Mb'] = infoDf['Tumor_Sample_Barcode'].apply(lambda x:  nmut_MbDict[x] if x in nmut_MbDict else None)
infoDf['pid'] = infoDf['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
infoDf['cancer_type'] = infoDf['pid'].apply(lambda x: cDict[x] if x in cDict else None)

#CREATE THE COHORTS!
signatureThreshold = .2 #the threshold at which we call a signature as existing
notHighThresh = 15
highThresh = 50

#POLE is divided into two classes: POLE endometrial and other pole
#If a case is 'mixed' POLE plus MMR I consider it POLE
poleEndometrialIds = set(impactSigs[(impactSigs['cancer_type'] == 'Endometrial Cancer') & (impactSigs['Nmut_Mb'] > highThresh) 
& ((impactSigs['Signature.10'] > signatureThreshold) | (impactSigs['Signature.14'] > signatureThreshold))]['Tumor_Sample_Barcode'])
    
otherPoleIds = set(impactSigs[(impactSigs['cancer_type'] != 'Endometrial Cancer') & (impactSigs['Nmut_Mb'] > highThresh) 
& ((impactSigs['Signature.10'] > signatureThreshold) | (impactSigs['Signature.14'] > signatureThreshold))]['Tumor_Sample_Barcode'])

#MMR is divided into three classes: Colorectal MMR, Endometrial MMR and Other MMR
colorectalMMRIds = set(impactSigs[(impactSigs['cancer_type'] == 'Colorectal Cancer') & (impactSigs['Signature.MMR'] > signatureThreshold) 
& ((impactSigs['Signature.10'] < signatureThreshold) & (impactSigs['Signature.14'] < signatureThreshold))]['Tumor_Sample_Barcode'])

endometrialMMRIds = set(impactSigs[(impactSigs['cancer_type'] == 'Endometrial Cancer') & (impactSigs['Signature.MMR'] > signatureThreshold) 
& ((impactSigs['Signature.10'] < signatureThreshold) & (impactSigs['Signature.14'] < signatureThreshold))]['Tumor_Sample_Barcode'])
 
otherMMRIds = set(impactSigs[(impactSigs['cancer_type'] != 'Colorectal Cancer') & (impactSigs['cancer_type'] != 'Endometrial Cancer') & (impactSigs['Nmut_Mb'] > highThresh)  & (impactSigs['Signature.MMR'] > signatureThreshold)
& ((impactSigs['Signature.10'] < signatureThreshold) | (impactSigs['Signature.14'] < signatureThreshold))]['Tumor_Sample_Barcode'])

#TMZ hypermutated is a specific group
gliomaTMZIds = set(impactSigs[(impactSigs['cancer_type'] == 'Glioma') & (impactSigs['Nmut_Mb'] > highThresh) 
& (impactSigs['Signature.11'] > signatureThreshold)]['Tumor_Sample_Barcode'])

#other High groups are just based on TMB and cancer type
cancerTypesForControls = set(['Colorectal Cancer', 'Non-Small Cell Lung Cancer', 'Glioma', 'Melanoma', 'Endometrial Cancer', 'Bladder Cancer'])    

infoDf['cohort'] = infoDf.apply(lambda row:
    
    row['cancer_type'] + '__not_high' if row['Nmut_Mb'] < notHighThresh and row['cancer_type'] in cancerTypesForControls
    else 'POLE_Endometrial' if row['Tumor_Sample_Barcode'] in poleEndometrialIds
    else 'POLE_Other' if row['Tumor_Sample_Barcode'] in otherPoleIds
    else 'MMR_Other' if row['Tumor_Sample_Barcode'] in otherMMRIds
    else 'MMR_Colorectal' if row['Tumor_Sample_Barcode'] in colorectalMMRIds
    else 'MMR_Endometrial' if row['Tumor_Sample_Barcode'] in endometrialMMRIds
    else 'MMR_Other' if row['Tumor_Sample_Barcode'] in otherMMRIds
    else 'TMZ_glioma' if row['Tumor_Sample_Barcode'] in gliomaTMZIds
    else row['cancer_type'] + '_high' if row['Nmut_Mb'] > highThresh and row['cancer_type'] in set(['Non-Small Cell Lung Cancer', 'Melanoma', 'Bladder Cancer'])
    else None
    ,axis=1)
 
#OLD WAY    
#highMutBurdenThresh = 20
#hypermutatorMutBurdenThresh = 80
#infoDf['cohort'] = infoDf.apply(lambda row: 
#    row['cancer_type'] + '__not_high' if row['Nmut_Mb'] < 10
#    else row['cancer_type'] + '__high' if row['Nmut_Mb'] < 80
#    else row['cancer_type'] + '__hypermutant' if row['Nmut_Mb'] >= 80
#    else row['cancer_type'] + '__not_high' #simplify because almost all cases without an nmut_mb estimate are not hypermutated #TODO FIX
#    , axis = 1)
#orderingDict1 = {'Endometrial Cancer': 0, 'Colorectal Cancer': 1, 'Glioma': 2, 'Melanoma': 3, 'Bladder Cancer': 4, 'Non-Small Cell Lung Cancer':5}
#orderingDict2 = {'hypermutant': .1, 'high': .2, 'not_high': .3}
#infoDf['orderingVal'] = infoDf['cohort'].apply(lambda x: orderingDict1[x.split('__')[0]] + orderingDict2[x.split('__')[1]]) 

orderingValDict = {'POLE_Other': 0, 'MMR_Other': 1,
                   'POLE_Endometrial':2, 'MMR_Endometrial':3, 'Endometrial Cancer__not_high':4,
                   'MMR_Colorectal':5, 'Colorectal Cancer__not_high':6,
                   'TMZ_glioma':7, 'Glioma__not_high':8,
                   'Melanoma_high':9, 'Melanoma__not_high':10,
                   'Bladder Cancer_high':11, 'Bladder Cancer__not_high':12,
                   'Non-Small Cell Lung Cancer_high':13, 'Non-Small Cell Lung Cancer__not_high':14,
                   }

infoDf['orderingVal'] = infoDf['cohort'].apply(lambda x: orderingValDict[x] if x in orderingValDict else None)
    

#NO LONGER NEED TO DO THIS COHORTS HAVE THIS BAKED IN
###ADD DOMINANT SIGNATURE INFORMATION
#sigNames = [i for i in impactSigs.columns.values if 'Signature.' in i]
#impactSigs['dominantSignature'] = impactSigs.apply(lambda row: mutationSigUtils.get_dominant_signature(row.to_dict(), cols=sigNames, prefix='Signature.'), axis=1)
#dominantSigDict = dict(zip(impactSigs['Tumor_Sample_Barcode'], impactSigs['dominantSignature']))
#infoDf['dominantSignature'] = infoDf['Tumor_Sample_Barcode'].apply(lambda x: dominantSigDict[x] if x in dominantSigDict else None)
#infoDf['dominantSignature'] = infoDf.apply(lambda row: 'Not Enough Mutations' if row['Nmut'] < 10 else row['dominantSignature'], axis=1)  
#infoDf['signatureAetiology'] = infoDf['dominantSignature'].apply(lambda x:
#    'Smoking' if x == 'Signature.SMOKING'
#    else 'MMR' if x == 'Signature.MMR'
#    else 'APOBEC' if x == 'Signature.APOBEC'
#    else 'Age' if x == 'Signature.1'
#    else 'POLE' if x == 'Signature.10'
#    else 'Mixed POLE/MMR' if x == 'Signature.14'
#    else 'UV' if x == 'Signature.7'
#    else 'BRCA' if x == 'Signature.3'
#    else 'TMZ' if x == 'Signature.11'
#    else x if x == 'Not Enough Mutations'
#    else 'Other Signature'
#    )

#Age of diagnosis stuff
ageAtSequencing = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/clinical_msk_impact_repository_data_adj.txt')
ageAtSequencingDict = dict(zip(ageAtSequencing['#Sample Identifier'], ageAtSequencing['Age at Which Sequencing was Reported (Days)']))
    
infoDf['ageAtSequencing'] = infoDf['Tumor_Sample_Barcode'].apply(lambda x: ageAtSequencingDict[x] if x in ageAtSequencingDict else None)

###ADD more information about doublets
#summaryDfDoubleMuts = double_mutation_analysis_util.create_double_mutation_summary_vanilla_maf(mafLimitedToCancerType)
#nGeneDoubleWithOncogenicDict = dict(zip(summaryDfDoubleMuts['Tumor_Sample_Barcode'], summaryDfDoubleMuts['nGenesDoubleOncogenicPerCase']))
#infoDf['nGenesWithDoubleOncogenic'] = infoDf['Tumor_Sample_Barcode'].apply(lambda x: nGeneDoubleWithOncogenicDict[x] if x in nGeneDoubleWithOncogenicDict else None) 

infoDf['cancer_type_fill'] = infoDf.apply(lambda row: 'Other' if row['cohort'] == 'MMR_Other' or row['cohort'] == 'POLE_Other' else row['cancer_type'], axis=1)
infoDf[infoDf['cohort'].notnull()].to_csv('~/Desktop/WORK/dataForLocalPlotting/mutburdenBoxplotV2.tsv', sep='\t', index=False)



##################SIGNATURE COHORT INFO


impactSigs = pd.read_table(pathPrefix + '/ifs/res/taylorlab/impact_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
impactSigs['pid'] = impactSigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
renameDict = {key:value for (key,value) in [('mean_' + str(i), 'Signature.' + str(i)) for i in range(1,31)]}
impactSigs = impactSigs.rename(columns=renameDict)
impactSigs = mutationSigUtils.merge_signature_columns(impactSigs, mode='Stratton', drop=False, smokingMerge=True, confidence=False, mean=True, prefix='Signature.')

sigNames = [i for i in impactSigs.columns.values if 'Signature.' in i]
impactSigs['dominantSignature'] = impactSigs.apply(lambda row: mutationSigUtils.get_dominant_signature(row.to_dict(), cols=sigNames, prefix='Signature.'), axis=1)
impactSigs['cancer_type'] = impactSigs['pid'].apply(lambda x: cDict[x] if x in cDict else None)
dominantSigDict = dict(zip(impactSigs['Tumor_Sample_Barcode'], impactSigs['dominantSignature']))
impactSigs['dominantSignature'] = impactSigs.apply(lambda row: 'Not Enough Mutations' if row['Nmut'] < 10 else row['dominantSignature'], axis=1)  
impactSigs['signatureAetiology'] = impactSigs['dominantSignature'].apply(lambda x:
    'Smoking' if x == 'Signature.SMOKING'
    else 'MMR' if x == 'Signature.MMR'
    else 'APOBEC' if x == 'Signature.APOBEC'
    else 'Age' if x == 'Signature.1'
    else 'POLE' if x == 'Signature.10'
    else 'Mixed POLE/MMR' if x == 'Signature.14'
    else 'UV' if x == 'Signature.7'
    else 'BRCA' if x == 'Signature.3'
    else 'TMZ' if x == 'Signature.11'
    else x if x == 'Not Enough Mutations'
    else 'Other Signature'
    )
    
impactSigs.to_csv('~/Desktop/WORK/dataForLocalPlotting/signaturePlotting.tsv', index=False, sep='\t')



###
###
#########
##############
####################
#WORK FOR HOTSPOT FREQUENCY ETC

#ranks hotspots by their prevalence per gene
def assign_hotspot_ranking_dict(df):
    d = {}
    for index, row in df.iterrows():
        localD = {}
        refAminoAcid = row['ref']
        gene = row['Hugo_Symbol']
        position = row['Amino_Acid_Position']
        for entry in row['Var_AA'].split('|'):
            fullAltName = ''
            altAminoAcid, count = entry.split(':')
            localD[altAminoAcid] = count
        print localD
        return
    return d

def assign_hotspot_freq_dict(df):
    d = {}
    for index, row in df.iterrows():
        
        refAminoAcid = row['ref']
        gene = row['Hugo_Symbol']
        position = row['Amino_Acid_Position']
        for entry in row['Var_AA'].split('|'):
            fullAltName = ''
            altAminoAcid, count = entry.split(':')
            fullAltName = gene + ':' + refAminoAcid + position + altAminoAcid
            d[fullAltName] = float(count)/47000
    return d

#TODO MAKE SURE SPLICE HOTSPOTS ARE PROPERLY LABELED
def asses_hotspot_freqs(focusMaf, refMaf, incidenceD):
    recurrentTumorSupressors, recurrentOncogenes = maf_analysis_utils.enumerate_recurrently_mutated_tumor_supressors_and_oncogenes(refMaf, thresh=.05)
    allTumorSuppressors = analysis_utils.get_tumor_supressor_genes()
    
    hotspots = focusMaf[focusMaf['is-a-hotspot'] == 'Y']
    hotspots['geneAlt'] = hotspots.apply(lambda row: row['Hugo_Symbol'] + ':'+ row['HGVSp_Short'].strip('p.'), axis=1)
    hotspots['incidence'] = hotspots['geneAlt'].apply(lambda x: incidenceD[x] if x in incidenceD else None)
    hotspots['geneClass'] = hotspots['Hugo_Symbol'].apply(lambda x:
        'RecurrentTumorSupressor' if x in recurrentTumorSupressors
        else 'RecurrentOncogene' if x in recurrentOncogenes
        else 'Tumor Supressor' if x in allTumorSuppressors
        else 'Oncogene')
    hotspots['xOrderingVal'] = hotspots['geneClass'].apply(lambda x: 1 if x == 'RecurrentOncogene'
            else 2 if x == 'Oncogene'
            else 3 if x == 'RecurrentTumorSupressor'
            else 4)
    return hotspots

#OTHER MODE DOING ANALYSIS OF HOTSPOT FREQUENCIES

mafWithInfo = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/annotatedOncoPlusHotspotMafAllImpact_trinuc')
impactSigs = pd.read_table(pathPrefix + '/ifs/res/taylorlab/impact_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
impactSigs['pid'] = impactSigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/mskImpactAsOfMarch2019/dmp/mskimpact/data_clinical_sample.txt')
impactSigs['cancer_type'] = impactSigs['pid'].apply(lambda x: cDict[x] if x in cDict else None)
mafWithInfo['pid'] = mafWithInfo['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
mafWithInfo['cancer_type'] = mafWithInfo['pid'].apply(lambda x: cDict[x] if x in cDict else None)

renameMapping = {'Pleural Mesothelioma, Epithelioid Type': 'Mesothelioma',
                 'Breast Invasive Ductal Carcinoma': 'Breast Cancer',
                 'Bladder Urothelial Carcinoma': 'Bladder Cancer',
                 'Upper Tract Urothelial Carcinoma': 'Bladder Cancer',
                 'Colon Adenocarcinoma': 'Colorectal Cancer',
                 'Glioblastoma Multiforme': 'Glioma',
                 'Adenocarcinoma of the Gastroesophageal Junction': 'Esophagogastric Cancer',
                 'Pancreatic Neuroendocrine Tumor': 'Pancreatic Cancer',
                 'Endometrial Carcinoma': 'Endometrial Cancer',
                 'Stomach Adenocarcinoma': 'Esophagogastric Cancer',
                 'Rectal Adenocarcinoma': 'Colorectal Cancer',
                 'High-Grade Serous Ovarian Cancer': 'Ovarian Cancer',
                 'Breast Invasive Lobular Carcinoma': 'Breast Cancer',
                 'Oligodendroglioma': 'Glioma',
                 'Serous Ovarian Cancer': 'Ovarian Cancer',
                 'Prostate Adenocarcinoma': 'Prostate Cancer',
                 'Breast Invasive Carcinoma, NOS': 'Breast Cancer',
                 'Esophageal Adenocarcinoma': 'Esophagogastric Cancer',
                 'Invasive Breast Carcinoma': 'Breast Cancer',
                 'Pancreatic Adenocarcinoma': 'Pancreatic Cancer',
                 'Uterine Endometrioid Carcinoma': 'Endometrial Cancer',
                 'Colorectal Adenocarcinoma': 'Colorectal Cancer',
                 'Mucinous Adenocarcinoma of the Colon and Rectum': 'Colorectal Cancer'
                 }

impactSigs['cancer_type'] = impactSigs['cancer_type'].apply(lambda x: renameMapping[x] if x in renameMapping else x)

impactSigs.to_csv('~/Desktop/WORK/dataForLocalPlotting/sigsWithCType.tsv', index=False, sep='\t')








hotspotsDf = pd.read_table(pathPrefix + '/home/gavrilae/snp_output_final_pancan.txt')

hypermutationThresh = 50
endometrialHyperIds = set(impactSigs[(impactSigs['cancer_type'] == 'Endometrial Cancer') & (impactSigs['Nmut_Mb'] > hypermutationThresh)]['Tumor_Sample_Barcode'])
endometrialNotHyperIds = set(impactSigs[(impactSigs['cancer_type'] == 'Endometrial Cancer') & (impactSigs['Nmut_Mb'] <= hypermutationThresh)]['Tumor_Sample_Barcode'])
colorectalHyperIds = set(impactSigs[(impactSigs['cancer_type'] == 'Colorectal Cancer') & (impactSigs['Nmut_Mb'] > hypermutationThresh)]['Tumor_Sample_Barcode'])
colorectalNotHyperIds = set(impactSigs[(impactSigs['cancer_type'] == 'Colorectal Cancer') & (impactSigs['Nmut_Mb'] <= hypermutationThresh)]['Tumor_Sample_Barcode'])
gliomaHyperIds = set(impactSigs[(impactSigs['cancer_type'] == 'Glioma') & (impactSigs['Nmut_Mb'] > hypermutationThresh)]['Tumor_Sample_Barcode'])
gliomaNotHyperIds = set(impactSigs[(impactSigs['cancer_type'] == 'Glioma') & (impactSigs['Nmut_Mb'] <= hypermutationThresh)]['Tumor_Sample_Barcode'])


hotspotIncidenceD = assign_hotspot_freq_dict(hotspotsDf)

assign_hotspot_ranking_dict(hotspotsDf)

hotspotsInfoDfEndometrial = asses_hotspot_freqs(mafWithInfo[mafWithInfo['Tumor_Sample_Barcode'].isin(endometrialHyperIds)], mafWithInfo[mafWithInfo['Tumor_Sample_Barcode'].isin(endometrialNotHyperIds)], hotspotIncidenceD)
hotspotsInfoDfColorectal = asses_hotspot_freqs(mafWithInfo[mafWithInfo['Tumor_Sample_Barcode'].isin(colorectalHyperIds)], mafWithInfo[mafWithInfo['Tumor_Sample_Barcode'].isin(colorectalNotHyperIds)], hotspotIncidenceD)
hotspotsInfoDfGlioma = asses_hotspot_freqs(mafWithInfo[mafWithInfo['Tumor_Sample_Barcode'].isin(gliomaHyperIds)], mafWithInfo[mafWithInfo['Tumor_Sample_Barcode'].isin(gliomaNotHyperIds)], hotspotIncidenceD)

hotspotsInfoDfEndoReg = asses_hotspot_freqs(mafWithInfo[mafWithInfo['Tumor_Sample_Barcode'].isin(endometrialNotHyperIds)], mafWithInfo[mafWithInfo['Tumor_Sample_Barcode'].isin(endometrialNotHyperIds)], hotspotIncidenceD)
hotspotsInfoDfColorectalReg = asses_hotspot_freqs(mafWithInfo[mafWithInfo['Tumor_Sample_Barcode'].isin(colorectalNotHyperIds)], mafWithInfo[mafWithInfo['Tumor_Sample_Barcode'].isin(colorectalNotHyperIds)], hotspotIncidenceD)

hotspotsInfoDfEndometrial.to_csv('~/Desktop/WORK/dataForLocalPlotting/Endometrial_hotspotPrevalenceInfo.tsv', index=False, sep='\t')
#hotspotsInfoDfEndometrialRef.to_csv('~/Desktop/WORK/dataForLocalPlotting/hotspotPrevalenceInfo.tsv', index=False, sep='\t')

hotspotsInfoDfColorectal.to_csv('~/Desktop/WORK/dataForLocalPlotting/Colorectal_hotspotPrevalenceInfo.tsv', index=False, sep='\t')
#hotspotsInfoDfColorectalReg.to_csv('~/Desktop/WORK/dataForLocalPlotting/hotspotPrevalenceInfo.tsv', index=False, sep='\t')
hotspotsInfoDfGlioma.to_csv('~/Desktop/WORK/dataForLocalPlotting/Glioma_hotspotPrevalenceInfo.tsv', index=False, sep='\t')



