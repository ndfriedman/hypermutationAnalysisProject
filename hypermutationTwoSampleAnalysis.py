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

def get_samples_with_multiple_cases(df):
    if 'pid' not in df.columns.values:
        df['pid'] = df['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
    multipleSampleIds = set()
    df = df.drop_duplicates(subset=['Tumor_Sample_Barcode'])
    
    multipleSampleIds = set(df[df.groupby('pid')['pid'].transform('size') > 1]['pid'])
    
    return multipleSampleIds

#classifies mutations as occuring before, after or in both samples
def classify_muts():
    return 0

#a limited function that classifies all mutations as being private/shared in series data
def classify_muts_private_shared_only(maf):
    listOfDfs = []
    for pid in set(maf['pid']):
        pidMaf = maf[maf['pid'] == pid]
        pidMaf['duplicateCol'] = pidMaf.duplicated(subset=['variantUUID'], keep=False)
        pidMaf['sharedMut'] = pidMaf['duplicateCol'].apply(lambda x: 'shared' if x == True else 'private')
        listOfDfs.append(pidMaf)
    return pd.concat(listOfDfs)
  

mutationSigs = pd.read_table(pathPrefix + '/ifs/res/taylorlab/impact_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
mutationSigs['pid'] = mutationSigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])

colorectalDf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/subsettedMafs/Colorectal_HypermutantCaseMuts_MAF_ANNO.maf')
endometrialDf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/subsettedMafs/Endometrial_HypermutantCaseMuts_MAF_ANNO.maf')
gliomaDf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/subsettedMafs/Glioma_HypermutantCaseMuts_MAF_ANNO.maf')

colorectalMultipleSampleHypermutantPids = get_samples_with_multiple_cases(colorectalDf)
endometrialMultipleSampleHypermutantPids = get_samples_with_multiple_cases(endometrialDf)
gliomaMultipleSampleHypermutantPids = get_samples_with_multiple_cases(gliomaDf)








#########OLD VERSION
      
mutationSigs = pd.read_table(pathPrefix + '/ifs/res/taylorlab/impact_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
mutationSigs['pid'] = mutationSigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
mutationSigs['cancer_type'] = mutationSigs['pid'].apply(lambda x: cDict[x] if x in cDict else None)
hypermutantTids = set(mutationSigs[(mutationSigs['Nmut_Mb'] > 50) & (mutationSigs['cancer_type'].isin(set(['Glioma', 'Colorectal Cancer', 'Endometrial Cancer'])))]['Tumor_Sample_Barcode'])
hypermutantPids = set(mutationSigs[(mutationSigs['Nmut_Mb'] > 50) & (mutationSigs['cancer_type'].isin(set(['Glioma', 'Colorectal Cancer', 'Endometrial Cancer'])))]['pid'])
#hypermutantTids = set(mutationSigs[mutationSigs['Nmut_Mb'] > 100]['Tumor_Sample_Barcode'])
#hypermutantPids = set(mutationSigs[mutationSigs['Nmut_Mb'] > 100]['pid'])
multipleSampleHypermutantPids = get_samples_with_multiple_cases(mutationSigs[mutationSigs['Tumor_Sample_Barcode'].isin(hypermutantTids)])



mafWithClonalityInfo['pid'] = mafWithClonalityInfo['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
mafWithClonalityInfoMultiple = mafWithClonalityInfo[mafWithClonalityInfo['pid'].isin(multipleSampleHypermutantPids)]
mafWithClonalityInfoMultiple['variantUUID'] = mafWithClonalityInfoMultiple.apply(lambda row: str(row['Chromosome']) + '_' + str(row['Start_Position']) + '_' + str(row['Reference_Allele']) + '->' + str(row['Tumor_Seq_Allele2']), axis=1)           

mafWithClonalityInfoMultipleOnco = mafWithClonalityInfoMultiple[mafWithClonalityInfoMultiple['oncogenic'].notnull()]

#a lot of the way this data frame is set up is lazy, be wary for errors

df = classify_muts_private_shared_only(mafWithClonalityInfoMultiple)
sharedM, sharedLower, sharedUpper = analysis_utils.mean_confidence_interval(df[df['sharedMut'] == 'shared']['clonal'], confidence=0.95)
privateM, privateLower, privateUpper = analysis_utils.mean_confidence_interval(df[df['sharedMut'] == 'private']['clonal'], confidence=0.95)

#Add some specific information about oncogenecity as well
oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
df['oncogenicVariant'] = df['oncogenic'].apply(lambda x: True if x in oncogenicMutColNames else False)
oncoMaf = df[df['oncogenicVariant'] == True]
sharedMOnc, sharedLowerOnc, sharedUpperOnc = analysis_utils.mean_confidence_interval(oncoMaf[oncoMaf['sharedMut'] == 'shared']['clonal'], confidence=0.95)
privateMOnc, privateLowerOnc, privateUpperOnc = analysis_utils.mean_confidence_interval(oncoMaf[oncoMaf['sharedMut'] == 'private']['clonal'], confidence=0.95)

#manually make a dataframe to put this information in logically
dfToWrite = pd.DataFrame([
        {'class': 'private', 'mean': privateM, 'lowerBound': privateLower, 'upperBound': privateUpper},
        {'class': 'shared', 'mean': sharedM, 'lowerBound': sharedLower, 'upperBound': sharedUpper},
        {'class': 'private_oncogenic', 'mean': privateMOnc, 'lowerBound': privateLowerOnc, 'upperBound': privateUpperOnc},
        {'class': 'shared_oncogenic', 'mean': sharedMOnc, 'lowerBound': sharedLowerOnc, 'upperBound': sharedUpperOnc}
        ])

#has to be a csv cause of some dataframe error
dfToWrite.to_csv('~/Desktop/WORK/dataForLocalPlotting/multipleSampleClonalityPlotting.csv', index=False, sep=',')


mutIds = set(mutationSigs[(mutationSigs['cancer_type'] == 'Glioma') & (mutationSigs['Nmut_Mb'] > 50)]['Tumor_Sample_Barcode'])
muts = mafWithClonalityInfo[mafWithClonalityInfo['Tumor_Sample_Barcode'].isin(mutIds)]

print np.nanmean(muts['ccf_Mcopies_upper'])



mafWithClonalityInfoMultipleOnco['ccf_1copy_upper']





