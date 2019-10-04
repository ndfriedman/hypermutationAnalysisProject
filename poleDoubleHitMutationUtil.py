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
import itertools


sigs = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/hypermutationAnalysisProjectData/hypermutationAnalysisProjectSignatures.tsv')
sigs['pid'] = sigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
sigs['cancer_type'] = sigs['pid'].apply(lambda x: cDict[x] if x in cDict else None)
cancerTypeDetailedDict = analysis_utils.get_cancer_type_detailed_information(cancerTypeDfPath = pathPrefix + '/ifs/work/taylorlab/friedman/dmp/mskimpact/data_clinical_sample.txt', mode='pid')
sigs['cancer_type_detailed'] = sigs['pid'].apply(lambda x: cancerTypeDetailedDict[x] if x in cancerTypeDetailedDict else None)
detailedCancerTypeDict = dict(zip(sigs['Tumor_Sample_Barcode'], sigs['cancer_type_detailed']))

sigNames = [i for i in sigs.columns.values if 'mean_' in i]
sigs['dominantSignature'] = sigs.apply(lambda row: mutationSigUtils.get_dominant_signature(row.to_dict(), cols=sigNames, prefix='mean_'), axis=1)
dominantSigDict = dict(zip(sigs['Tumor_Sample_Barcode'], sigs['dominantSignature']))

maf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/hypermutationAnalysisProjectData/hypermutationAnalysisProjectMaf.tsv')
maf['vaf'] = maf.apply(lambda row: 1.0*row['t_alt_count']/(row['t_alt_count'] + row['t_ref_count']), axis = 1)

endometrialHypermutatorIds = set(sigs[(sigs['cancer_type'] == 'Endometrial Cancer') & (sigs['Nmut_Mb'] > 80)]['Tumor_Sample_Barcode'])
endometrialNotHypermutatorIds = set(sigs[(sigs['cancer_type'] == 'Endometrial Cancer') & (sigs['Nmut_Mb'] < 80)]['Tumor_Sample_Barcode'])

endometrialHypermutatorMuts = maf[maf['Tumor_Sample_Barcode'].isin(endometrialHypermutatorIds)]
endometrialHypermutatorMuts['quadNuc'] = endometrialHypermutatorMuts.apply(lambda row: mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)
endometrialHypermutatorMuts['dominantSignature'] = endometrialHypermutatorMuts['Tumor_Sample_Barcode'].apply(lambda x: dominantSigDict[x] if x in dominantSigDict else None)

mafAnnoMaf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/All.dmp_somatic_data_mutations_unfiltered.mafAnno.oncokb.hotspots.maf', sep='\t')
mafAnnoMaf['vaf'] = mafAnnoMaf.apply(lambda row: 1.0*row['t_alt_count']/(row['t_alt_count'] + row['t_ref_count']) if row['t_alt_count'] + row['t_ref_count'] > 0 else None, axis = 1)

endometrialHyperMafAnno = mafAnnoMaf[mafAnnoMaf['Tumor_Sample_Barcode'].isin(endometrialHypermutatorIds)]

fgaData = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/alexPFgaCopy.txt')
fgaDict = dict(zip(fgaData['Tumor_Sample_Barcode'], fgaData['FGA']))
mafAnnoMaf['fga'] = mafAnnoMaf['Tumor_Sample_Barcode'].apply(lambda x: fgaDict[x] if x in fDict else None)

endometrialLowFGANotHyper = mafAnnoMaf[mafAnnoMaf['Tumor_Sample_Barcode'].isin(endometrialNotHypermutatorIds) & (mafAnnoMaf['fga'] < .1) & (mafAnnoMaf['fga'].notnull())]

mafToAnalyze = endometrialLowFGANotHyper

mafToAnalyze = endometrialHyperMafAnno 

for case in set(mafToAnalyze['Tumor_Sample_Barcode']):
    caseMaf = mafToAnalyze[mafToAnalyze['Tumor_Sample_Barcode'] == case]
    
    stdVaf = np.std(caseMaf['vaf'])
    maxVaf = max(caseMaf['vaf'])
    meanVaf = np.nanmean(caseMaf['vaf'])
    medianVaf = np.nanmedian(caseMaf['vaf'])
    veryHighAfMaf = caseMaf[caseMaf['vaf'] > 2* (max(meanVaf, medianVaf))]
    
    #print case 
    #print 'nMuts: ', caseMaf.shape[0]
    #print 'max vaf: ', maxVaf
    #print 'mean vaf: ', meanVaf
    #print 'median vaf', medianVaf
    #print 'standard deviation: ', stdVaf
    #print 'n mutations above 2*max(mean, median): ', veryHighAfMaf.shape[0]
    
    genesAtVeryHiVaf = set(veryHighAfMaf['Hugo_Symbol'])
    for gene in genesAtVeryHiVaf:
        geneMaf = caseMaf[caseMaf['Hugo_Symbol'] == gene]
        startPos = list(geneMaf['Start_Position'])
        pairs = [(x,y) for i,x in enumerate(startPos) for j,y in enumerate(startPos) if i != j]
        closeToEachOther = [x for x in pairs if abs(x[0] - x[1]) < 100]
        if len(closeToEachOther) > 1:
            positions = set(list(itertools.chain(*closeToEachOther)))
            positionMaf = geneMaf[geneMaf['Start_Position'].isin(positions)]
            
            print case
            print gene, meanVaf, maxVaf
            print positionMaf['vaf']
            print positionMaf['Start_Position']
            print ''
            
            #print case, gene, closeToEachOther
#abs(x - y) < 1000
    
print case, sigs[sigs['Tumor_Sample_Barcode'] == 'P-0020378-T01-IM6']['mean_25']
    
    

for case in set(endometrialHypermutatorMuts['Tumor_Sample_Barcode']):
    caseMaf = endometrialHypermutatorMuts[endometrialHypermutatorMuts['Tumor_Sample_Barcode'] == case]
    poleMuts = caseMaf[caseMaf['Hugo_Symbol'] == 'POLE']
    if poleMuts.shape[0] > 0:
        maxPoleVaf = max(poleMuts['vaf'])
        mutsDoubleOverPole = caseMaf[caseMaf['vaf'] >= 2*maxPoleVaf]
        if mutsDoubleOverPole.shape[0] < 15 and mutsDoubleOverPole.shape[0] > 0:
            print case
            print list(mutsDoubleOverPole['Hugo_Symbol'])
            print set(mutsDoubleOverPole['quadNuc'])
            print list(mutsDoubleOverPole['dominantSignature'])[0]
    else:
        pass
    

    
    
    