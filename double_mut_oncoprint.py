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
import double_mutation_analysis_util


def assign_gene_status(row, gene):
    if row[gene +'_lohPlusOncogenic'] == True and row[gene +'_doubleOncogenic'] == True: return 'Double Oncogenic & LOH' 
    elif row[gene +'_lohPlusOncogenic'] == True: return 'Oncogenic+LOH' 
    elif row[gene +'_doubleOncogenic'] == True : return 'DoubleOncogenic'
    elif row[gene +'_doubleWithOncogenic'] == True: return 'Oncogenic+Other'
    elif row[gene +'_oncogenicMut'] == True: return 'OncogenicMut'
    elif row[gene +'_loh'] == True: return 'LohOnly' 
    else: return 'Neither oncogenic mut nor loh'
        

sigs = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/hypermutationAnalysisProjectData/hypermutationAnalysisProjectSignatures.tsv')
sigs['pid'] = sigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
sigs['cancer_type'] = sigs['pid'].apply(lambda x: cDict[x] if x in cDict else None)
cancerTypeDetailedDict = analysis_utils.get_cancer_type_detailed_information(cancerTypeDfPath = pathPrefix + '/ifs/work/taylorlab/friedman/dmp/mskimpact/data_clinical_sample.txt', mode='pid')
sigs['cancer_type_detailed'] = sigs['pid'].apply(lambda x: cancerTypeDetailedDict[x] if x in cancerTypeDetailedDict else None)

#TWO DIFFERENT WAYS TO DO THIS
#ALL but....
endometrialHypermutatorSubtypeIds = set(sigs[(sigs['cancer_type'] == 'Endometrial Cancer')
 & (sigs['cancer_type_detailed'] != 'Uterine Serous Carcinoma/Uterine Papillary Serous Carcinoma')
 & (sigs['cancer_type_detailed'] != 'Uterine Carcinosarcoma/Uterine Malignant Mixed Mullerian Tumor')]['Tumor_Sample_Barcode'])
#ENdometrioid Only
endometrialHypermutatorSubtypeIds = set(sigs[(sigs['cancer_type'] == 'Endometrial Cancer') & 
 ((sigs['cancer_type_detailed'] == 'Uterine Endometrioid Carcinoma') | (sigs['cancer_type_detailed'] == 'Endometrial Carcinoma') | (sigs['cancer_type_detailed'] == 'Uterine Mixed Endometrial Carcinoma'))
    ]['Tumor_Sample_Barcode'])

endometrialIds = set(sigs[sigs['cancer_type'] == 'Endometrial Cancer']['Tumor_Sample_Barcode'])
    
mutBurdenDict = dict(zip(sigs['Tumor_Sample_Barcode'], sigs['Nmut_Mb']))

mafAnnoMaf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/All.dmp_somatic_data_mutations_unfiltered.mafAnno.oncokb.hotspots.maf', sep='\t')
mafAnnoMaf['isHypermutator'] = mafAnnoMaf['Tumor_Sample_Barcode'].apply(lambda x: True if x in hypermutatorIds else False)

endometrialHypermutatorPossibleSubtype = mafAnnoMaf[mafAnnoMaf['Tumor_Sample_Barcode'].isin(endometrialHypermutatorSubtypeIds)]
endometrialMaf = mafAnnoMaf[mafAnnoMaf['Tumor_Sample_Barcode'].isin(endometrialIds)]

summaryDfDoubleMuts = double_mutation_analysis_util.create_double_mutation_summary_vanilla_maf(endometrialHypermutatorPossibleSubtype, genes = set(['TP53', 'ARID1A', 'PIK3CA', 'PTEN', 'CTNNB1', 'NF1']))
summaryDfLOH = maf_analysis_utils.do_gene_loh_summary(endometrialHypermutatorPossibleSubtype, genes = set(['TP53', 'ARID1A', 'PIK3CA', 'PTEN', 'CTNNB1', 'NF1']))

mergedDf = pd.merge(summaryDfDoubleMuts, summaryDfLOH, on='Tumor_Sample_Barcode')
mergedDfLOHInfoPresent = mergedDf[mergedDf['Tumor_Sample_Barcode'].isin(set(endometrialHypermutatorPossibleSubtype['Tumor_Sample_Barcode']))]
mergedDfLOHInfoPresent['Nmut_Mb'] = mergedDfLOHInfoPresent['Tumor_Sample_Barcode'].apply(lambda x: mutBurdenDict[x] if x in mutBurdenDict else None)

mergedDfLOHInfoPresent['PTEN_status'] = mergedDfLOHInfoPresent.apply(lambda row: assign_gene_status(row, 'PTEN'), axis=1)
mergedDfLOHInfoPresent['PIK3CA_status'] = mergedDfLOHInfoPresent.apply(lambda row: assign_gene_status(row, 'PIK3CA'), axis=1)
mergedDfLOHInfoPresent['ARID1A_status'] = mergedDfLOHInfoPresent.apply(lambda row: assign_gene_status(row, 'ARID1A'), axis=1)
mergedDfLOHInfoPresent['TP53_status'] = mergedDfLOHInfoPresent.apply(lambda row: assign_gene_status(row, 'TP53'), axis=1)
mergedDfLOHInfoPresent['NF1_status'] = mergedDfLOHInfoPresent.apply(lambda row: assign_gene_status(row, 'NF1'), axis=1)
mergedDfLOHInfoPresent['CTNNB1_status'] = mergedDfLOHInfoPresent.apply(lambda row: assign_gene_status(row, 'CTNNB1'), axis=1)

#ADD SIGNATURE INFO
sigs['dominantSignature'] = sigs.apply(lambda row: mutationSigUtils.get_dominant_signature(row.to_dict()), axis=1)
dominantSigDict = dict(zip(sigs['Tumor_Sample_Barcode'], sigs['dominantSignature']))

discardSignatureNmutMbThresh = 8 #samples with 8 nmut_mb have between 9 and 13 impact mutations
mergedDfLOHInfoPresent['dominantSignature'] = mergedDfLOHInfoPresent['Tumor_Sample_Barcode'].apply(lambda x: dominantSigDict[x] if x in dominantSigDict else None)
mergedDfLOHInfoPresent['signatureAetiology'] = mergedDfLOHInfoPresent.apply(lambda row:
    'Not enough mutations' if row['Nmut_Mb'] < discardSignatureNmutMbThresh
    else 'POLE'  if row['dominantSignature'] == 'mean_10'
    else 'AGE'  if row['dominantSignature'] == 'mean_1'
    else 'MIXED POLE MMR' if row['dominantSignature'] == 'mean_14'
    else 'MMR' if (row['dominantSignature'] == 'mean_6') | (row['dominantSignature'] == 'mean_20') | (row['dominantSignature'] == 'mean_26') | (row['dominantSignature'] == 'mean_15')
    else 'Other'
    ,axis=1)


mergedDfLOHInfoPresent.to_csv('~/Desktop/WORK/dataForLocalPlotting/endometrialOncoprintData.tsv', sep='\t', index=False)



mergedDfLOHInfoPresent.columns.values

