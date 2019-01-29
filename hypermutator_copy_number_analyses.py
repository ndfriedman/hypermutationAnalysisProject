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

fgaData = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/alexPFgaCopy.txt')
fgaDict = dict(zip(fgaData['Tumor_Sample_Barcode'], fgaData['FGA']))

sigs = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/hypermutationAnalysisProjectData/hypermutationAnalysisProjectSignatures.tsv')
sigs['FGA'] = sigs['Tumor_Sample_Barcode'].apply(lambda x: fgaDict[x] if x in fgaDict else None)
sigs['pid'] = sigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
sigs['cancer_type'] = sigs['pid'].apply(lambda x: cDict[x] if x in cDict else None)
sigs = mutationSigUtils.merge_signature_columns(sigs)

sigsWithFGAData = sigs[sigs['FGA'].notnull()]

print np.nanmean(sigsWithFGAData[sigsWithFGAData['Nmut_Mb'] > 200]['FGA'])

print np.nanmean(sigsWithFGAData[sigsWithFGAData['cancer_type'] == 'Breast Cancer']['FGA'])

print np.nanmean(sigsWithFGAData[(sigsWithFGAData['cancer_type'] == 'Endometrial Cancer') & (sigsWithFGAData['Nmut_Mb'] > 140) & (sigsWithFGAData['Nmut_Mb'] < 1400)]['FGA'])

print np.nanmean(sigsWithFGAData[(sigsWithFGAData['cancer_type'] == 'Endometrial Cancer') & (sigsWithFGAData['mean_14'] > .25)]['FGA'])

endometrial = sigsWithFGAData[sigsWithFGAData['cancer_type'] == 'Endometrial Cancer']
endometrial.to_csv('~/Desktop/WORK/dataForLocalPlotting/endometrialCNAData.tsv', sep='\t', index=False)

bladder = sigsWithFGAData[sigsWithFGAData['cancer_type'] == 'Bladder Cancer']
bladder.to_csv('~/Desktop/WORK/dataForLocalPlotting/bladderCNAData.tsv', sep='\t', index=False)

colorectal = sigsWithFGAData[sigsWithFGAData['cancer_type'] == 'Colorectal Cancer']
colorectal.to_csv('~/Desktop/WORK/dataForLocalPlotting/colorectalCNAData.tsv', sep='\t', index=False)

glioma = sigsWithFGAData[sigsWithFGAData['cancer_type'] == 'Glioma']
glioma.to_csv('~/Desktop/WORK/dataForLocalPlotting/gliomaCNAData.tsv', sep='\t', index=False)

melanoma = sigsWithFGAData[sigsWithFGAData['cancer_type'] == 'Melanoma']
melanoma.to_csv('~/Desktop/WORK/dataForLocalPlotting/melanomaCNAData.tsv', sep='\t', index=False)



