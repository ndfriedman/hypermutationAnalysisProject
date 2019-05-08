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
import mutation_modeling_util


#CALCULATES THE DIFFERENTIAL EXPECTED oncogenic mut burden as if each signature were dominant alone
def calculate_oncogenic_mut_susceptibility_of_genes_by_signature(oncogenicSDict):
    listOfDicts = []
    sigNames = ['Signature.' + str(i) for i in range(1,31)]
    for i in range(1,31):
        curSig = 'Signature.' + str(i)
        d = {}
        for s in sigNames:
            d[s] = 0
        d[curSig] = 1
        #PRETEND we got a case with 100% signature i on the decomposition
        quadNucFractions = mutation_modeling_util.get_quadnuc_fracs_given_decomposition(d, spectraPath = pathPrefix + '/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt')
        v = mutation_modeling_util.get_expected_oncogenic_val_given_quadnuc_fractions(quadNucFractions, oncogenicSDict, 'IMPACT_468')
        listOfDicts.append({'Signature_Name': curSig, 'ExpectedFracOfMutsOncogenic': v})
    return pd.DataFrame(listOfDicts)

def expand_data_for_plot(infoDict, n=1250):
    listOfDicts = []
    for i in range(1,n):
        if i%50==0:print i
        nmut_mbIM6 = (i*1000000.0)/1139322
        for key, value in infoDict.items():
            listOfDicts.append({'Signature': key, 'Nmut_Expected': i*value, 'Nmut_Mb': nmut_mbIM6})
    return pd.DataFrame(listOfDicts)

allPossibleImpactMuts = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/panImpactSimulatedMutationDataSummary.tsv')
oncogenicSusceptibilityDict = mutation_modeling_util.calculate_quadnuc_based_oncogenic_susceptibility_dict(allPossibleImpactMuts)
geneSpecificOncogenicSusceptibilityDict = mutation_modeling_util.calculate_quadnuc_and_gene_based_oncogenic_susceptibility_dict(allPossibleImpactMuts)

#get the oncogenic fraction by signature and get MMR/indel fraction info
oncogenicFracDf = calculate_oncogenic_mut_susceptibility_of_genes_by_signature(oncogenicSusceptibilityDict)
mmrIndelOncogenicity = mutation_modeling_util.calculate_pan_impact_likelihood_of_oncogenic_mmr_indel(repeatRegionInfo = pathPrefix + '/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/IMPACT_repeat_stats.txt')
mmr6OncSusceptibility = float(oncogenicFracDf[oncogenicFracDf['Signature_Name'] == 'Signature.6']['ExpectedFracOfMutsOncogenic'].iloc[0])
#we model mmr oncogenic susceptibility by assuming 1/3 of TMB comes from indels and 2/3 come from snps
mmrOncSusceptibility = (1.0/3)*mmrIndelOncogenicity + (2.0/3)*mmr6OncSusceptibility

#add this to our data frame (oncogenicFracDf)
oncogenicFracDf = oncogenicFracDf.append(pd.DataFrame([{'ExpectedFracOfMutsOncogenic': mmrOncSusceptibility, 'Signature_Name': 'MSI_(MMR+1/3indels)'}]))
oncogenicFracDict = dict(zip(oncogenicFracDf['Signature_Name'], oncogenicFracDf['ExpectedFracOfMutsOncogenic']))
df = expand_data_for_plot(oncogenicFracDict)
df['colorOpaque'] = df['Signature'].apply(lambda x: 1 if x in set(['Signature.28', 'Signature.10', 'MSI_(MMR+1/3indels)']) else .2)

df.to_csv('~/Desktop/WORK/dataForLocalPlotting/expectedOncogenicDataBySig.tsv', index=False, sep='\t')


#TODO--plot the mapped opberved number of oncogenic mutations by dominant signature across bands

##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################


allImpactFilteredMuts = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/tempScriptFiles/filteredMafAnnotated.maf')
allImpactFilteredMuts = maf_analysis_utils.fix_mll_genes(allImpactFilteredMuts)
oncogenicMutBurdenSummaryDict = maf_analysis_utils.get_n_oncogenic_muts_per_case_dict(allImpactFilteredMuts)

nmutMbInfo = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/tmbInfo.tsv')
nmutMbDict = dict(zip(nmutMbInfo['Tumor_Sample_Barcode'], nmutMbInfo['Nmut_Mb']))

#$#$#$
sigsInfo = pd.read_table(pathPrefix + '/ifs/res/taylorlab/impact_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
sigsInfo['pid'] = sigsInfo['Tumor_Sample_Barcode'].apply(lambda x: x[:9])

sigsInfo = mutationSigUtils.merge_signature_columns(sigsInfo, mode='Stratton', drop=False, smokingMerge=True, confidence=False, mean=True, prefix='Signature.')
sigsInfo['dominantSig'] = sigsInfo.apply(lambda row: mutationSigUtils.get_dominant_signature(row.to_dict(), cols=None, prefix='mean'), axis=1)
sigsToInclude = set(['mean_1', 'mean_2', 'mean_3', 'mean_4', 'mean_6', 'mean_7', 'mean_10', 'mean_11', 'mean_14', 'mean_15', 'mean_17', 'mean_20', 'mean_21', 'mean_26', 'insufficientMutBurden'])
sigsInfo['dominantSigAdj'] = sigsInfo['dominantSig'].apply(lambda x: x if x in sigsToInclude else 'Other')

sigsInfo['oncogenicMutsPerCase'] = sigsInfo['Tumor_Sample_Barcode'].apply(lambda x: oncogenicMutBurdenSummaryDict[x] if x in oncogenicMutBurdenSummaryDict else None)
sigsInfo = sigsInfo[sigsInfo['oncogenicMutsPerCase'].notnull()]

sigsInfo['Nmut_Mb'] = sigsInfo['Tumor_Sample_Barcode'].apply(lambda x: nmutMbDict[x] if x in nmutMbDict else None)
sigsInfo = sigsInfo[sigsInfo['Nmut_Mb'].notnull()]

mmrSigs = set(['mean_6', 'mean_15', 'mean_20', 'mean_21', 'mean_26'])
sigsInfo['sigDisplayLabel'] = sigsInfo['dominantSigAdj'].apply(lambda x: 'mmrRelated' if x in mmrSigs
        else 'AGE' if x == 'mean_1' else 'APOBEC' if x in set(['mean_2', 'mean_13']) else 'POLE' if x == 'mean_10' else 'BRCA' if x == 'mean_3' else 'UV' if x == 'mean_7'
        else 'POLE_MMR' if x == 'mean_14' else 'TMZ' if x == 'mean_11' else 'SMOKING' if x == 'mean_4' else 'other')

sigsInfo.to_csv('~/Desktop/WORK/dataForLocalPlotting/observedOncogenicDataBySig.tsv', index=False, sep='\t')



