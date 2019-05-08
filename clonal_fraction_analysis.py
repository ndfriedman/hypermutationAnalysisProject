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
import clonality_analysis_util

def adjust_dict(data, minNCases=10):
    newDict = {}
    otherVals = []
    for key, value in data.items():
        if len(value) >= minNCases:
            newDict[key] = value #add error bars
        else:
            otherVals = otherVals + value
    newDict['otherGenes'] = otherVals
    return newDict

"""def summarize_vaf_info(data, minNCases=10):
  d = {}
   otherMutDiffList = []
    for key, value in data.items():
        if len(value) >= minNCases:
            d[key] = np.nanmean(value) #add error bars
        else:
            otherMutDiffList = otherMutDiffList + value
    d['otherGenes'] = np.nanmedian(otherMutDiffList)
    return d"""

def fix_mll_genes(maf):
    maf['Hugo_Symbol'] = maf['Hugo_Symbol'].apply(lambda x:
        'KMT2A' if x == 'MLL'
        else 'KMT2B' if x == 'MLL2'
        else 'KMT2C' if x == 'MLL3'
        else x)   
    return maf

def process_data(df, n=20):
    #FILTERING on oncogenic mutations with an annotated protein change
    df = df[df['HGVSp_Short'].notnull()]
    df = df[df['oncogenic'].notnull()]
    
    d = clonality_analysis_util.calculate_delta_vaf_across_mutation_pairs(df)
    dDeltaVaf = dict()
    dMaxVaf = dict()
    for key, value in d.items():
        dDeltaVaf[key] = value[0]
        dMaxVaf[key] = value[1]
     
    dDeltaDictAdj = adjust_dict(dDeltaVaf, minNCases=10)
    dMaxVafAdj = adjust_dict(dMaxVaf, minNCases=10)
    
    listOfDicts = []
    for key in dDeltaDictAdj.keys():
        m1, dVafLower, dVafUpper = analysis_utils.mean_confidence_interval(dDeltaDictAdj[key], confidence=0.95)
        m2, maxVafLower, maxVafUpper = analysis_utils.mean_confidence_interval(dMaxVafAdj[key], confidence=0.95)
        
        listOfDicts.append({'gene': key,
                            'dVafLower': dVafLower, 'dVafUpper': dVafUpper,
                            'dVafMedian': np.nanmedian(dDeltaDictAdj[key]), 'dVafMean': m1,
                            'maxVafLower': maxVafLower, 'maxVafUpper': maxVafUpper,
                            'maxVafMedian': np.nanmedian(dMaxVafAdj[key]), 'maxVafMean': m2
                            })
    
    df = pd.DataFrame(listOfDicts)
    return df
        

colorectalDf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/subsettedMafs/Colorectal_HypermutantCaseMuts_MAF_ANNO.maf')
colorectalDf = fix_mll_genes(colorectalDf)
dfColorectal = process_data(colorectalDf)
dfColorectal.to_csv('~/Desktop/WORK/dataForLocalPlotting/colorectalVafInfo.tsv', index=False, sep='\t')


endometrialDf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/subsettedMafs/Endometrial_HypermutantCaseMuts_MAF_ANNO.maf')
endometrialDf = fix_mll_genes(endometrialDf)
dfEndometrial = process_data(endometrialDf, n = 30)
dfEndometrial.to_csv('~/Desktop/WORK/dataForLocalPlotting/endometrialVafInfo.tsv', index=False, sep='\t')

#sharedM, sharedLower, sharedUpper = analysis_utils.mean_confidence_interval(df[df['sharedMut'] == 'shared']['clonal'], confidence=0.95)
"""dColorectal = clonality_analysis_util.calculate_delta_vaf_across_mutation_pairs(colorectalDf[colorectalDf['oncogenic'].notnull()])
#averagesColorectal = [(key, np.nanmean(value)) for key, value in dColorectal.items() if len(value) > 7]
averagesColorectal = summarize_delta_vaf_info(dColorectal, minNCases=10)
sortedL = sorted(averagesColorectal, key = lambda x: x[1])
"""

####SET UP THE FILES to run MAF-anno on FOR THIS ANALYSIS
df = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/tempScriptFiles/attempt2Maf.maf')
colorectalHyperMaf = df[df['Tumor_Sample_Barcode'].isin(analysis_utils.get_ids_by_hypermutant_status(hypermutantIdDir=pathPrefix + '/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/hypermutationStatusIds', cancerType='Colorectal Cancer', hypermutantStatus = 'Hypermutated'))]
endometrialHyperMaf = df[df['Tumor_Sample_Barcode'].isin(analysis_utils.get_ids_by_hypermutant_status(hypermutantIdDir=pathPrefix + '/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/hypermutationStatusIds', cancerType='Endometrial Cancer', hypermutantStatus = 'Hypermutated'))]
gliomaHyperMaf = df[df['Tumor_Sample_Barcode'].isin(analysis_utils.get_ids_by_hypermutant_status(hypermutantIdDir=pathPrefix + '/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/hypermutationStatusIds', cancerType='Glioma', hypermutantStatus = 'Hypermutated'))]

colorectalHyperMaf.to_csv(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/subsettedMafs/Colorectal_HypermutantCaseMuts.maf', index=False, sep='\t')
endometrialHyperMaf.to_csv(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/subsettedMafs/Endometrial_HypermutantCaseMuts.maf', index=False, sep='\t')
gliomaHyperMaf.to_csv(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/subsettedMafs/Glioma_HypermutantCaseMuts.maf', index=False, sep='\t')


print len(set(gliomaHyperMaf['Tumor_Sample_Barcode']))





