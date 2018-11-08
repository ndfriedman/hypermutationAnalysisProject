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



mafWithHotspotInfo = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/annotatedOncoPlusHotspotMafAllImpact_trinuc')


Counter(mafWithHotspotInfo[mafWithHotspotInfo['is-a-hotspot'] == 'Y']['oncogenic'])

#summaryInfoPanCohort = maf_analysis_utils.summarize_per_case_mutation_info_for_mafs(mafWithHotspotInfo)
summaryInfoPanCohort.to_csv(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/mafSummaryData/mafOncogenicMutSummaryData.tsv', sep='\t', index=False)
#summaryInfoPanCohort = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/mafSummaryData/mafOncogenicMutSummaryData.tsv', sep='\t', index=False)

summaryInfoPanCohort.columns.values
summaryInfoPanCohort['activatingMutToNmutRatio'] = summaryInfoPanCohort.apply(lambda row: 1.0*row['nActivatingMutations']/row['Nmut'],axis = 1)
summaryInfoPanCohort['hotspotMutToNmutRatio'] =  summaryInfoPanCohort.apply(lambda row: 1.0*row['nHotspotMuts']/row['Nmut'],axis = 1)
summaryInfoPanCohort['oncogenicMutToNmutRatio'] =  summaryInfoPanCohort.apply(lambda row: 1.0*row['nOncogenicMuts']/row['Nmut'],axis = 1)

np.nanmean(summaryInfoPanCohort[(summaryInfoPanCohort['Nmut'] < 10000) & (summaryInfoPanCohort['Nmut'] > 1000)]['activatingMutToNmutRatio'])






