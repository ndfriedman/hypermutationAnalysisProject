#written by Noah Friedman (a template for scripts to be excuted in the spyder environment
#A script to analayze clonality of oncogenic mutations in hypermutator cases
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

def compare_clonality_across_mut_groups(
        mutsWithClonalInfo, #a maf with mutations with clonal info
        caseIds, #ids for the cases that fit our class
        orderedGenesImplicatedInDisease, #an ordered list of the 15 most recurrently mutated genes in a cancer type
        sigMotifDict): # a dictionary mapping each of the 96 quadnucs to a label
    
    caseMuts = mutsWithClonalInfo[mutsWithClonalInfo['Tumor_Sample_Barcode'].isin(caseIds)]
    top5Genes = set(orderedGenesImplicatedInDisease[:5])
    fiveThrough15Genes = set(orderedGenesImplicatedInDisease[5:])
    
    
    caseMuts['geneClass'] = caseMuts['Hugo_Symbol'].apply(lambda x:
        ('-').join(orderedGenesImplicatedInDisease[:5]) if x in top5Genes
        else ('-').join(orderedGenesImplicatedInDisease[5:]) if x in fiveThrough15Genes
        else 'other genes')
    
    caseMuts['mutEffectLabel'] = caseMuts.apply(lambda row:
        'hotspot' if row['is-a-hotspot'] == 'Y'
        else 'oncogenic but not hotspot' if (str(row['oncogenic']) == 'Likely Oncogenic' or str(row['oncogenic']) == 'Oncogenic' or str(row['oncogenic']) == 'Predicted Oncogenic')
        else 'no annotation'
        , axis=1)
    
    caseMuts['sigMotifLabel'] = caseMuts['quadNuc'].apply(lambda x:
        sigMotifDict[x] if x != None else None)
        
    
    return caseMuts


        
   
def init_trinuc_mapping_dict(trinucGroups):
    d = {}
    for entry in trinucGroups:
        trinucs, label = entry
        for trinuc in trinucs:
            d[trinuc] = label
    return d    

mafWithOncogenicInfo = mafHere
#ALEX g's big maf
mafWithClonalityInfo = pd.read_table(pathPrefix + '/ifs/res/taylorlab/ang46/ext/dmp/mskimpact/mutation_data.txt')
mafWithClonalityInfoClonalOnly = mafWithClonalityInfo[mafWithClonalityInfo['clonal'].notnull()]
mafWithClonalityInfoClonalOnly['quadNuc'] = mafWithClonalityInfoClonalOnly.apply(lambda row: mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)
#mafWithClonalityInfo = mafBackground

impactSigs = pd.read_table(pathPrefix + '/ifs/res/taylorlab/impact_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
impactSigs['pid'] = impactSigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
impactSigs['cancer_type'] = impactSigs['pid'].apply(lambda x: cDict[x] if x in cDict else None)
allTrinucs = analysis_utils.get_all_possible_quadNucs()

ids = impactSigs[(impactSigs['Nmut'] > 100) & (impactSigs['cancer_type'] != 'Colorectal Cancer') & (impactSigs['cancer_type'] != 'Endometrial Cancer') & (impactSigs['mean_10'] > .25)]['cancer_type']
for v in ids:
    print v

impactSigs[impactSigs['Tumor_Sample_Barcode'] == 'P-0005270-T02-IM6']['mean_6']

###################************************#######################
#GLIOMA analyses
gliomaSigs = impactSigs[impactSigs['cancer_type'] == 'Glioma']
#gliomasWithTMZIds = set(gliomaSigs[(gliomaSigs['mean_11'] >= .25) & (gliomaSigs['Nmut'] >= 20)]['Tumor_Sample_Barcode'])
gliomasWithTMZIds = set(gliomaSigs[gliomaSigs['Nmut'] <= 20]['Tumor_Sample_Barcode'])

orderedGenesImplicatedInGlioma = ['TERT', 'TP53', 'IDH1', 'PTEN', 'EGFR', 'ATRX', 'NF1', 'CIC', 'PIK3CA', 'PIK3CR1', 'RB1', 'NOTCH1', 'PDGFRA', 'BRAF', 'PTPN11']
tmzTop = set(['ACTC', 'CCTC', 'GCTC', 'TCTC'])
tmzSecond = set(['ACTT', 'CCTT', 'GCTT', 'TCTT'])
tmzThird = set(['ACTA', 'CCTA', 'GCTA', 'TCTA'])
tmzFourth = set(['ACTG', 'CCTG', 'GCTG', 'TCTG'])
notTMZ = allTrinucs - tmzTop - tmzSecond - tmzThird - tmzFourth
trinucMappings = [(tmzTop, 'topTMZ'),
                  (tmzSecond, 'secondTMZ'),
                  (tmzThird, 'thirdTMZ'),
                  (tmzFourth, 'fourthTMZ'),
                  (notTMZ, 'other')]
tmzTrinucsDict = init_trinuc_mapping_dict(trinucMappings)
gliomaClonalDf = compare_clonality_across_mut_groups(mafWithClonalityInfoClonalOnly, gliomasWithTMZIds, orderedGenesImplicatedInGlioma, tmzTrinucsDict)

gliomaClonalDfAdj = gliomaClonalDf[['geneClass', 'mutEffectLabel', 'sigMotifLabel', 'Hugo_Symbol', 'clonal']] #make it little because of some errors in the overall dataframe in R
#ALERT Idk if we want this
gliomaClonalDfAdj = gliomaClonalDfAdj[gliomaClonalDfAdj['sigMotifLabel'].notnull()]

gliomaClonalDfAdj.to_csv('~/Desktop/dataForLocalPlotting/gliomaClonalityData.tsv', sep='\t', index=False)


np.nanmean(gliomaClonalDfAdj[gliomaClonalDfAdj['geneClass'] != 'TERT-TP53-IDH1-PTEN-EGFR']['clonal'])

######################################################################
#POLE Endometrial Analyses

#PROBLEM!!!!!!!!!! MOST ENDOMETRIAL DATA DOES NOT HAVE CLONALITY ESTIMATES

endometrialSigs = impactSigs[impactSigs['cancer_type'] == 'Endometrial Cancer']
endometrialWithPOLEIds = set(endometrialSigs[(endometrialSigs['mean_10'] >= .25) & (endometrialSigs['Nmut'] >= 20)]['Tumor_Sample_Barcode'])

endoDf = mafWithClonalityInfo[mafWithClonalityInfo['Tumor_Sample_Barcode'].isin(endometrialWithPOLEIds)]
Counter(endoDf[endoDf['clonal'].isnull()])

len(set(endoDf['Tumor_Sample_Barcode']))
#compare_clonality_across_mut_groups(mafWithClonalityInfoClonalOnly, gliomasWithTMZIds, orderedGenesImplicatedInGlioma, tmzTrinucsDict)

####################################################################





#COLORECTAL HYPER Analyses
#THESE cases are mostly MMR with a smattering of POLE and friends

colorectalSigs = impactSigs[impactSigs['cancer_type'] == 'Colorectal Cancer']
#ALERT currently mmr 6 only
colorectalHyper = colorectalSigs[(colorectalSigs['Nmut'] > 100) & (colorectalSigs['mean_6'] > .25)]
colorectalHyperIds = set(colorectalHyper['Tumor_Sample_Barcode'])
colorectalHyperMuts = mafWithClonalityInfoClonalOnly[mafWithClonalityInfoClonalOnly['Tumor_Sample_Barcode'].isin(colorectalHyperIds)]
orderedGenesImplicatedInColorectal = ['APC', 'TP53', 'KRAS', 'PIK3CA', 'SMAD4', 'FBXW7', 'TCF7L2', 'BRAF', 'SOX9', 'KMT2D', 'ARID1A', 'PTPRS', 'RNF43', 'ATM', 'PTPRT']
mmrMotifs = set(['CCAT', 'ACTG', 'CCTG', 'GCTA', 'GCTC', 'GCTG', 'GCTT'])
trinucMappings = [(mmrMotifs, 'MMR'),
        (allTrinucs - mmrMotifs, 'other')]
colorectalTrinucsDict = init_trinuc_mapping_dict(trinucMappings)
colorectalClonalDf = compare_clonality_across_mut_groups(mafWithClonalityInfoClonalOnly, colorectalHyperIds, orderedGenesImplicatedInColorectal, colorectalTrinucsDict)
colorectalClonalDfAdj = colorectalClonalDf[['geneClass', 'mutEffectLabel', 'sigMotifLabel', 'Hugo_Symbol', 'clonal']] #make it little because of some errors in the overall dataframe in R
#TEMP
colorectalClonalDfAdj = colorectalClonalDfAdj[colorectalClonalDfAdj['sigMotifLabel'].notnull()]
colorectalClonalDfAdj.to_csv('~/Desktop/dataForLocalPlotting/colorectalClonalityData.tsv', sep='\t', index=False)

np.nanmean(colorectalClonalDf[colorectalClonalDf['sigMotifLabel'] == 'MMR']['clonal'])
######################################################################






#MELANOMA Analyses
melanomaSigs = impactSigs[impactSigs['cancer_type'] == 'Melanoma']
melanomaHyperUV = melanomaSigs[(melanomaSigs['Nmut'] > 100) & (melanomaSigs['mean_7'] > .25)]
melanomaHyperUVIds = set(melanomaHyperUV['Tumor_Sample_Barcode'])
melanomaHyperMuts = mafWithClonalityInfoClonalOnly[mafWithClonalityInfoClonalOnly['Tumor_Sample_Barcode'].isin(melanomaHyperUVIds)]
orderedGenesImplicatedInMelanoma = ['TERT', 'BRAF', 'PTPRT', 'NRAS', 'NF1', 'GRIN2A', 'TP53', 'PTPRD', 'ROS1', 'PAK5', 'KMT2D', 'ERBB4', 'PIK3C2G', 'TP63', 'CDKN2A']
uvMotifs = set(['CCTA', 'CCTC', 'CCTT', 'TCTA', 'TCTC', 'TCTG', 'TCTT'])
trinucMappings = [(uvMotifs, 'UV'),
        (allTrinucs - uvMotifs, 'other')]
melanomaTrinucsDict = init_trinuc_mapping_dict(trinucMappings)
melanomaClonalDf = compare_clonality_across_mut_groups(mafWithClonalityInfoClonalOnly, melanomaHyperUVIds, orderedGenesImplicatedInMelanoma, melanomaTrinucsDict)
melanomaClonalDfAdj = melanomaClonalDf[['geneClass', 'mutEffectLabel', 'sigMotifLabel', 'Hugo_Symbol', 'clonal']] #make it little because of some errors in the overall dataframe in R
#TEMP
melanomaClonalDfAdj = melanomaClonalDfAdj[melanomaClonalDfAdj['sigMotifLabel'].notnull()]
melanomaClonalDfAdj.to_csv('~/Desktop/dataForLocalPlotting/melanomaClonalityData.tsv', sep='\t', index=False)

###################################################################






#LUNG analyses?
lungSigs = impactSigs[impactSigs['cancer_type'] == 'Non-Small Cell Lung Cancer']
lungHyperSmoking = lungSigs[(lungSigs['Nmut'] > 100) & (lungSigs['mean_4'] > .25)]
lungHyperIds = set(lungHyperSmoking['Tumor_Sample_Barcode'])
lungHyperMuts = mafWithClonalityInfoClonalOnly[mafWithClonalityInfoClonalOnly['Tumor_Sample_Barcode'].isin(lungHyperIds)]
orderedGenesImplicatedInLungCancer = ['TP53', 'KRAS', 'EGFR', 'STK11', 'KEAP1', 'RBM10', 'PTPRD', 'SMARCA4', 'KMT2D', 'NF1', 'FAT1', 'CDKN2A', 'ARID1A', 'PTPRT', 'ATM', 'ALK']
lungMotifs = set(['ACAA', 'TCAA', 'CCAA', 'GCAA', 
                  'ACAT', 'TCAT', 'CCAT', 'GCAT',
                  'ACAC', 'TCAC', 'CCAC', 'GCAC',
                  'ACAG', 'TCAG', 'CCAG', 'GCAG'])
trinucMappings = [(lungMotifs, 'Smoking'),
                  (allTrinucs- lungMotifs, 'other motif')]
lungCancerTrinucsDict = init_trinuc_mapping_dict(trinucMappings)
lungClonalDf = compare_clonality_across_mut_groups(mafWithClonalityInfoClonalOnly, lungHyperIds, orderedGenesImplicatedInLungCancer, lungCancerTrinucsDict)
lungClonalDfAdj = lungClonalDf[['geneClass', 'mutEffectLabel', 'sigMotifLabel', 'Hugo_Symbol', 'clonal']] #make it little because of some errors in the overall dataframe in R
#TEMP
lungClonalDfAdj = lungClonalDfAdj[lungClonalDfAdj['sigMotifLabel'].notnull()]

lungClonalDfAdj.to_csv('~/Desktop/dataForLocalPlotting/lungClonalityData.tsv', sep='\t', index=False)

#################################################################





#BLADDER analyses??
bladderSigs = impactSigs[impactSigs['cancer_type'] == 'Bladder Cancer']
bladderSigsHyper = bladderSigs[((bladderSigs['mean_2'] > .25) | (bladderSigs['mean_13'] > .25)) & (bladderSigs['Nmut'] > 100)]
bladderHyperIds = set(bladderSigsHyper['Tumor_Sample_Barcode'])
bladderHyperMuts = mafWithClonalityInfoClonalOnly[mafWithClonalityInfoClonalOnly['Tumor_Sample_Barcode'].isin(bladderHyperIds)]
orderedGenesImplicatedInBladderCancer = ['TERT', 'TP53', 'KDM6A', 'FGFR3', 'KMT2D', 'ARID1A', 'PIK3CA', 'CREBBP', 'CDKN1A', 'KMT2C', 'RB1', 'ERBB2', 'FAT1', 'EP300', 'ATM']
apobecMotifs = set(['TCGA', 'TCGC', 'TCGT', 'TCTA', 'TCTT', 'TCTC'])
trinucMappings = [(apobecMotifs, 'APOBEC'),
        (allTrinucs - apobecMotifs, 'other')]
bladderCancerTrinucsDict = init_trinuc_mapping_dict(trinucMappings)
bladderClonalDf = compare_clonality_across_mut_groups(mafWithClonalityInfoClonalOnly, bladderHyperIds, orderedGenesImplicatedInBladderCancer, bladderCancerTrinucsDict)
bladderClonalDfAdj = bladderClonalDf[['geneClass', 'mutEffectLabel', 'sigMotifLabel', 'Hugo_Symbol', 'clonal']] #make it little because of some errors in the overall dataframe in R
#TEMP
bladderClonalDfAdj = bladderClonalDfAdj[bladderClonalDfAdj['sigMotifLabel'].notnull()]

bladderClonalDfAdj.to_csv('~/Desktop/dataForLocalPlotting/bladderClonalityData.tsv', sep='\t', index=False)

np.nanmean(bladderClonalDfAdj[bladderClonalDfAdj['geneClass'] == 'TERT-TP53-KDM6A-FGFR3-KMT2D']['clonal'])

for case in set(colorectalHyperMuts['Tumor_Sample_Barcode']):
    caseMuts = colorectalHyperMuts[colorectalHyperMuts['Tumor_Sample_Barcode'] == case]
    print 1.0*caseMuts[caseMuts['clonal'] == 1.0].shape[0]/caseMuts.shape[0]



oncogenicDict = dict(zip(mafWithoncogenicInfoTMZGliomas['uniqueMutationIdentifier'], mafWithoncogenicInfoTMZGliomas['oncogenic']))

#TODO work with oncogenic variant annotations from the other MAF

#postTMZMutations = mafWithClonalityInfo[mafWithClonalityInfo['Tumor_Sample_Barcode'].isin(postTMZSampleIds)]
#postTMZMutations['uniqueMutationIdentifier'] = postTMZMutations.apply(lambda row: str(row['Hugo_Symbol'])+ '_' + str(row['HGVSp_Short']) + '_@' + str(row['Start_Position']), axis=1)

#clonalityInfoDict = dict(zip(postTMZMutations['uniqueMutationIdentifier'], postTMZMutations['clonal']))
#postTMZMutations['oncogenicBetter'] = postTMZMutations['uniqueMutationIdentifier'].apply(lambda x: oncogenicDict[x] if x in oncogenicDict else None)
#oncoKbOncogenicAnnotations = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
#mafWithoncogenicInfoTMZGliomas['clonal'] = mafWithoncogenicInfoTMZGliomas['uniqueMutationIdentifier'].apply(lambda x: clonalityInfoDict[x] if x in clonalityInfoDict else None)




    