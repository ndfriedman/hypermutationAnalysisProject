#written by Noah Friedman (a template for scripts to be excuted in the spyder environment
#CONTAINS DIFFERENT LITTLE UTILITIES TO DO HYPERMUTATN CLASSIFICATION
#ALSO RELIES ON plotAndDefineHypermutationThreshold.R
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
from scipy import stats

#Calculate the per case NMUT information
def calculate_per_case_nmut_info_filtered_maf(maf):
    listOfDicts = []
    cntr = 0
    for case in set(maf['Tumor_Sample_Barcode']):
        if cntr%100 ==0: print cntr
        caseMaf = maf[maf['Tumor_Sample_Barcode'] == case]
        nMuts = caseMaf.shape[0]
        listOfDicts.append({'Tumor_Sample_Barcode': case, 'Nmut': nMuts})
        cntr += 1
    return pd.DataFrame(listOfDicts)

#
def adjust_data_given_filtered_nmut_info(df, dictOfNmut, nmutFilteredThresh = 15):
    
    df['hypermutantClassification'] = df.apply(lambda row:
        row['hypermutantClassification'] if row['hypermutantClassification'] != 'Hypermutated'
        else 'indeterminateLowFilteredMutBurden' if row['Tumor_Sample_Barcode'] not in dictOfNmut
        else 'indeterminateLowFilteredMutBurden' if dictOfNmut[row['Tumor_Sample_Barcode']] <= nmutFilteredThresh
        else row['hypermutantClassification']
        , axis=1)
    return df

impactSigs = pd.read_table(pathPrefix + '/ifs/res/taylorlab/impact_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
impactSigs['pid'] = impactSigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
impactSigs['cancer_type'] = impactSigs['pid'].apply(lambda x: cDict[x] if x in cDict else None)
impactSigs.to_csv('~/Desktop/WORK/dataForLocalPlotting/sigsWithCType.tsv', index=False, sep='\t')

nmutFilteredDict = dict(zip(filteredMafNmutDf['Tumor_Sample_Barcode'], filteredMafNmutDf['Nmut']))

sigClassificationDir = pathPrefix + '/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/hypermutationStatusIds'
# a loop to classify information
"""for val in os.listdir(sigClassificationDir):
    print val,
    fullPath = os.path.join(sigClassificationDir, val)
    df = pd.read_table(fullPath)
    #df = adjust_data_given_filtered_nmut_info(df, nmutFilteredDict)
    #df.to_csv(fullPath, index=False, sep='\t')
"""   
#filteredMaf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/mskImpactAsOfMarch2019/dmp/mskimpact/data_mutations_extended.txt', skiprows=[0])
#filteredMafNmutDf = calculate_per_case_nmut_info_filtered_maf(filteredMaf) 
#filteredMafNmutDf.to_csv(pathPrefix + '/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/nmutFilteredInfo.tsv', index=False, sep='\t')

analysis_utils.get_ids_by_hypermutant_status(pathPrefix + '/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/hypermutationStatusIds', 
                                   'Endometrial Cancer', hypermutantStatus = 'Hypermutated')


