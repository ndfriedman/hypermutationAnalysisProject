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
import mutationSigUtils

def summarize_hypermutated_cases_and_dominant_sig(hypermutationStatusDir, dominantSigsDict):
    cancerTypesToExclude = set(['Bladder_Cancer', 'Cancer_of_Unknown_Primary', 'Non-Small_Cell_Lung_Cancer', 'Small_Cell_Lung_Cancer', 'Lung_Adenocarcinoma', 'Melanoma'])
    listOfDs = []
    for filename in os.listdir(hypermutationStatusDir):
        filepath = os.path.join(hypermutationStatusDir, filename)
        df = pd.read_table(filepath)
        hypermutants = df[df['hypermutantClassification'] == 'Hypermutated']
        cancerType = filename.strip('.tsv')
        if cancerType not in cancerTypesToExclude:
            print cancerType
            for index, row in hypermutants.iterrows():
                if row['Tumor_Sample_Barcode'] in dominantSigsDict:
                    listOfDs.append({'Tumor_Sample_Barcode': row['Tumor_Sample_Barcode'],
                               'Nmut_Mb': row['Nmut_Mb'],
                               'cancerType': cancerType,
                               'dominantSig': dominantSigsDict[row['Tumor_Sample_Barcode']]
                                })
    return pd.DataFrame(listOfDs)

sigsData = pd.read_table(pathPrefix + '/ifs/res/taylorlab/impact_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
sigsData = mutationSigUtils.merge_signature_columns(sigsData, mode='Stratton', drop=True, smokingMerge=False, confidence=True, mean=True, prefix='mean_')
sigsData['dominantSignature'] = sigsData.apply(lambda row: mutationSigUtils.get_dominant_signature(row.to_dict(), cols=None, prefix='mean'), axis=1)
dSigsDict = dict(zip(sigsData['Tumor_Sample_Barcode'], sigsData['dominantSignature']))

hypermutantStatusDir = pathPrefix + '/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/hypermutationStatusIds/'
df = summarize_hypermutated_cases_and_dominant_sig(hypermutantStatusDir, dSigsDict)
countsC = Counter(df['cancerType'])
df['orderingValCType'] = df['cancerType'].apply(lambda x: countsC[x])

df['dominantSig'] = df['dominantSig'].apply(lambda x: x if x in set(['mean_MMR', 'mean_1', 'mean_APOBEC', 'mean_10', 'mean_11', 'mean_14']) else 'other')
countsS = Counter(df['dominantSig'])
df['orderingValSig'] = df['dominantSig'].apply(lambda x: countsS[x] if x != 'other' else 1000)

df['cancerTypeAdj'] = df['cancerType'].apply(lambda x: x if x in set(['Endometrial_Cancer', 'Colorectal_Cancer', 'Glioma']) else 'other')

df.to_csv('~/Desktop/WORK/dataForLocalPlotting/hypermutantCohortInfo.tsv', sep='\t', index=False)



