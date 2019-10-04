#A FILE intended to contain assorted calls to utilities to make essential config files

#written by Noah Friedman (a template for scripts to be excuted in the spyder environment
import sys
import argparse
import os
import pandas as pd
import numpy as np
import re

from collections import Counter

pathPrefix = ''
if os.getcwd() == '/Users/friedman/Desktop/mnt':
	pathPrefix = '/Users/friedman/Desktop/mnt'

sys.path.append(pathPrefix + '/ifs/work/taylorlab/friedman/myUtils')
import analysis_utils 
import mutationSigUtils 
import maf_analysis_utils

configFileDir = pathPrefix + '/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles'

#i
#MAKE Indel fraction files
sigs = pd.read_table(pathPrefix + '/ifs/res/taylorlab/impact_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
sigs['dominantSignature'] = sigs.apply(lambda row: mutationSigUtils.get_dominant_signature(row.to_dict(), cols=None, prefix='mean', notEnoughMuts= True), axis=1)

impactFilteredMaf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/tempScriptFiles/filteredMafAnnotated_trinuc.maf')

minCasesToDeriveIndelFracEstimate = 25
listOfDicts = []
for sig in set(sigs['dominantSignature']):
    if sig != 'insufficientMutBurden':
        dominantSigCases = set(sigs[(sigs['dominantSignature'] == sig) &
                                    (sigs['Nmut_Mb'] > 10) #only take cases with a non-trivial mutation burden
                                    &(sigs[sig] > 1.0/3)] #only take cases where the dominant sig is greater than 33% of total signature burden
    ['Tumor_Sample_Barcode'])
        
        if len(dominantSigCases) > minCasesToDeriveIndelFracEstimate:
            dominantSigMaf = impactFilteredMaf[impactFilteredMaf['Tumor_Sample_Barcode'].isin(dominantSigCases)]
            listOfDicts.append({'Signature': re.sub('mean_', 'Signature.', sig), 
                                'Indel_Fraction': maf_analysis_utils.get_indel_frac_for_cohort(dominantSigMaf)})
    
panImpactIndelFrac = maf_analysis_utils.get_indel_frac_for_cohort(impactFilteredMaf)
listOfDicts.append({'Signature': 'All_IMPACT_Cases', 'Indel_Fraction': panImpactIndelFrac})
df = pd.DataFrame(listOfDicts)
df.to_csv(pathPrefix + '/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/indelRateBySignature.tsv', index=False, sep='\t')      
        
df.columns.values

#ii
#summarize msi area for indels
msiIndelInfo = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/IMPACT_repeat_stats.txt')
msiIndelInfo['msiBp'] = msiIndelInfo.apply(lambda row: row['bp']*row['rep.pct'], axis=1)
msiIndelInfo['msiProb'] = msiIndelInfo['msiBp'].apply(lambda x: 1.0*x/np.nansum(msiIndelInfo['bp']))
    
msiProbByGene = msiIndelInfo[['Hugo_Symbol', 'bp', 'msiProb']]
msiProbByGene.to_csv(os.path.join(configFileDir, 'relativeMsiProbabilityByGene.tsv'), index = False, sep='\t')





##########################################333











