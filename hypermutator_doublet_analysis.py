
#written by Noah Friedman (a template for scripts to be excuted in the spyder environment

import sys
import argparse
import os
import pandas as pd
import numpy as np
import scipy.stats

from collections import Counter




pathPrefix = ''
if os.getcwd() == '/Users/friedman/Desktop/mnt':
	pathPrefix = '/Users/friedman/Desktop/mnt'

sys.path.append(pathPrefix + '/ifs/work/taylorlab/friedman/myUtils')
import analysis_utils 
import mutationSigUtils 
import maf_analysis_utils

def make_pik3ca_summary_df(muts, ids):
    listOfDicts = []
    for case in ids:
        localD = {}
        caseMaf = muts[muts['Tumor_Sample_Barcode'] == case]
        pik3caMuts = caseMaf[caseMaf['Hugo_Symbol'] == 'PIK3CA']
        pik3caMutsNames = list(pik3caMuts['HGVSp_Short'])
        nPik3caMuts = len(pik3caMutsNames)
        pik3caHotspots = pik3caMuts[pik3caMuts['is-a-hotspot'] == 'Y']
        pik3caHotspotsNames =  list(pik3caHotspots['HGVSp_Short'])
        nPik3caHotspots = len(pik3caHotspotsNames)
        localD['nPik3ca'] = nPik3caMuts
        localD['pik3caNames'] = pik3caMutsNames
        localD['nPik3caHotspot'] = nPik3caHotspots
        localD['pik3caHotspotNames'] = pik3caHotspotsNames
        localD['Tumor_Sample_Barcode'] = case
        listOfDicts.append(localD)
    return pd.DataFrame(listOfDicts)

def summarize_n_double_mutations(mutsDf, adjNmutDict):
    ids = set(mutsDf['Tumor_Sample_Barcode'])
    listOfDicts = []
    cntr = 0
    for case in ids:
        if cntr%500 == 0: print cntr, len(ids)
        cntr += 1
        
        localD = {}
        caseMaf = mutsDf[mutsDf['Tumor_Sample_Barcode'] == case]
        hotspotMuts = caseMaf[caseMaf['is-a-hotspot'] == 'Y']
        hotspotGenes = list(hotspotMuts['Hugo_Symbol'])
        multipleGenesHotspots = [item for item, count in Counter(hotspotGenes).items() if count > 1]
        
        genesAll = list(caseMaf['Hugo_Symbol'])
        multipleGenesMuts = [item for item, count in Counter(genesAll).items() if count > 1]
        nMultipleGeneMuts = len(multipleGenesMuts)
        nmut = None
        if case in adjNmutDict:
            nmut = adjNmutDict[case]
        if nMultipleGeneMuts > 0 and nmut != None:
            ratioNmutToNGenesWithDoublets = 1.0*nmut/nMultipleGeneMuts
        else: 
            ratioNmutToNGenesWithDoublets = None
        
        localD['nGenesWithDoubletMutations'] = nMultipleGeneMuts
        localD['Nmut_Mb'] = nmut
        localD['nGenesWithDoubletHotspotMutations'] = len(multipleGenesHotspots)
        localD['genesWithDoubletHotspots'] = multipleGenesHotspots
        localD['nmutToNDoubletRatio'] = ratioNmutToNGenesWithDoublets
        localD['Tumor_Sample_Barcode'] = case
        localD['isHypermutator'] = caseMaf['isHypermutator'].iloc[0]
        
        listOfDicts.append(localD)
        
    return pd.DataFrame(listOfDicts)
        
        

# little utility function to do string logix
def contains_util(x, val):
    if val in str(x): return True
    else: return False
    

impactSigs = pd.read_table(pathPrefix + '/ifs/res/taylorlab/impact_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
impactSigs['pid'] = impactSigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
impactSigs['cancer_type'] = impactSigs['pid'].apply(lambda x: cDict[x] if x in cDict else None)
mafWithHotspotInfo = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/annotatedOncoPlusHotspotMafAllImpact_trinuc')
mafWithHotspotInfo['pid'] = mafWithHotspotInfo['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
mafWithHotspotInfo['cancer_type'] = mafWithHotspotInfo['pid'].apply(lambda x: cDict[x] if x in cDict else None)
#I found 8 mislabeled cancers fix those
poleEndometrialNotInMafIDs = set(['P-0027652-T01-IM6', 'P-0029690-T01-IM6', 'P-0030322-T01-IM6', 'P-0030372-T01-IM6', 'P-0032496-T01-IM6', 'P-0032548-T01-IM6', 'P-0033995-T01-IM6', 'P-0035147-T01-IM6'])
mafWithHotspotInfo['cancer_type'] = mafWithHotspotInfo['cancer_type'].apply(lambda x: 'Endometrial Cancer' if x in poleEndometrialNotInMafIDs else x)


#mafWithClonalityInfo = pd.read_table(pathPrefix + '/ifs/res/taylorlab/ang46/ext/dmp/mskimpact/mutation_data.txt')
#mafWithHotspotInfo = mafWithClonalityInfo


###############POLE ANALYSES
poleEndometrialIds = set(impactSigs[(impactSigs['cancer_type'] == 'Endometrial Cancer') & (impactSigs['mean_10'] > .25) & (impactSigs['Nmut'] > 100)]['Tumor_Sample_Barcode'])
poleEndometrialNotInMafIDs = set(['P-0027652-T01-IM6', 'P-0029690-T01-IM6', 'P-0030322-T01-IM6', 'P-0030372-T01-IM6', 'P-0032496-T01-IM6', 'P-0032548-T01-IM6', 'P-0033995-T01-IM6', 'P-0035147-T01-IM6'])
poleEndometrialIds = poleEndometrialIds | poleEndometrialNotInMafIDs

mmrEndoIds = set(impactSigs[(impactSigs['cancer_type'] == 'Endometrial Cancer') & (impactSigs['mean_6'] > .25) & (impactSigs['Nmut'] > 100) & (impactSigs['mean_10'] < .15)]['Tumor_Sample_Barcode'])

poleEndoMuts = mafWithHotspotInfo[mafWithHotspotInfo['Tumor_Sample_Barcode'].isin(poleEndometrialIds)]
mmrEndoMuts = mafWithHotspotInfo[mafWithHotspotInfo['Tumor_Sample_Barcode'].isin(mmrEndoIds)]
summarize_n_double_mutations(poleEndoMuts)
summarize_n_double_mutations(mmrEndoMuts)

poleEndometrialSummary = make_pik3ca_summary_df(poleEndoMuts, poleEndometrialIds)
set(poleEndoMuts[poleEndoMuts['clonal'] > 0]['Tumor_Sample_Barcode'])
endometrialNotHighMutBurdenIds = set(impactSigs[(impactSigs['Nmut'] <= 25) & (impactSigs['cancer_type'] == 'Endometrial Cancer')]['Tumor_Sample_Barcode'])
endometrialMuts = mafWithHotspotInfo[mafWithHotspotInfo['Tumor_Sample_Barcode'].isin(endometrialNotHighMutBurdenIds)]
endometrialSummary = make_pik3ca_summary_df(endometrialMuts, endometrialNotHighMutBurdenIds)

endometrialSummary = endometrialSummary[endometrialSummary['nPik3ca'] > 0]
poleEndometrialSummary = poleEndometrialSummary[poleEndometrialSummary['nPik3ca'] > 0]

#####################################
endometrialSummary['Tumor_Sample_Barcode']

print endometrialSummary[endometrialSummary['nPik3caHotspot'] > 1].shape
print poleEndometrialSummary[poleEndometrialSummary['nPik3caHotspot'] > 1].shape

print endometrialSummary[endometrialSummary['nPik3ca'] > 1].shape
print poleEndometrialSummary[poleEndometrialSummary['nPik3ca'] > 1].shape

print endometrialSummary[(endometrialSummary['nPik3ca'] > 1) & (endometrialSummary['nPik3caHotspot'] > 0)].shape
print poleEndometrialSummary[(poleEndometrialSummary['nPik3ca'] > 1) & (poleEndometrialSummary['nPik3caHotspot'] > 0)].shape

endometrialSummary['r88qPresent'] = endometrialSummary['pik3caHotspotNames'].apply(lambda x: contains_util(x, 'R88Q'))
poleEndometrialSummary['r88qPresent'] = poleEndometrialSummary['pik3caHotspotNames'].apply(lambda x: contains_util(x, 'R88Q'))

endometrialSummary['y1021cPresent'] = endometrialSummary['pik3caHotspotNames'].apply(lambda x: contains_util(x, 'Y1021C'))
poleEndometrialSummary['y1021cPresent'] = poleEndometrialSummary['pik3caHotspotNames'].apply(lambda x: contains_util(x, 'Y1021C'))

print endometrialSummary[endometrialSummary['r88qPresent'] == True].shape
print poleEndometrialSummary[poleEndometrialSummary['r88qPresent'] == True].shape

print endometrialSummary[endometrialSummary['y1021cPresent'] == True].shape
print poleEndometrialSummary[poleEndometrialSummary['y1021cPresent'] == True].shape

print 'N doublet hotspots'
print scipy.stats.fisher_exact([[21, 13], [216,13]], alternative='two-sided')

print 'N pik3ca doublets'
print scipy.stats.fisher_exact([[32, 22], [205,4]], alternative='two-sided')

print 'N pik3ca hotspot + second hit'
print scipy.stats.fisher_exact([[30, 20], [207,6]], alternative='two-sided')

print 'N r88q'
print scipy.stats.fisher_exact([[12, 10], [225,16]], alternative='two-sided')

print 'N y1021c'
print scipy.stats.fisher_exact([[0, 3], [227,23]], alternative='two-sided')

print impactSigs[(impactSigs['cancer_type'] == 'Breast Cancer') & (impactSigs['Nmut'] > 100)]['mean_APOBEC']
impactSigs[impactSigs['Nmut'] == 75]['Nmut_Mb']

print impactSigs[impactSigs['Tumor_Sample_Barcode'].isin(set(['P-0027652-T01-IM6', 'P-0029690-T01-IM6', 'P-0030322-T01-IM6', 'P-0030372-T01-IM6', 'P-0032496-T01-IM6', 'P-0032548-T01-IM6', 'P-0033995-T01-IM6', 'P-0035147-T01-IM6']))]['cancer_type']

#ADD a couple of cancer types and 
impactSigs['cancer_type'] = impactSigs['cancer_type'].apply(lambda x: 'Endometrial Cancer' if x in poleEndometrialNotInMafIDs else x)
##########################
N_MUT_MB_THRESHOLD = 40
impactSigs = mutationSigUtils.merge_signature_columns(impactSigs, mode='Stratton', drop=True)

poleHypermutatorIds = set(impactSigs[(impactSigs['mean_10'] > .25) & (impactSigs['Nmut_Mb'] > N_MUT_MB_THRESHOLD)]['Tumor_Sample_Barcode'])
mmrHypermutatorIds = set(impactSigs[(impactSigs['mean_MMR'] > .25) & (impactSigs['Nmut_Mb'] > N_MUT_MB_THRESHOLD) & (impactSigs['mean_10'] < .25)]['Tumor_Sample_Barcode']) #IGNORE POLE + MMR
apobecHypermutatorIds = set(impactSigs[(impactSigs['mean_APOBEC'] > .25) & (impactSigs['Nmut_Mb'] > N_MUT_MB_THRESHOLD)]['Tumor_Sample_Barcode'])
tmzHypermutatorIds = set(impactSigs[(impactSigs['mean_11'] > .25) & (impactSigs['Nmut_Mb'] > N_MUT_MB_THRESHOLD)]['Tumor_Sample_Barcode'])
uvHypermutatorIds = set(impactSigs[(impactSigs['mean_7'] > .25) & (impactSigs['Nmut_Mb'] > N_MUT_MB_THRESHOLD)]['Tumor_Sample_Barcode'])

#######################
sigsNmutDict = dict(zip(impactSigs['Tumor_Sample_Barcode'], impactSigs['Nmut_Mb']))

#############ENDOMETRIAL HYPER ANALYSES
"""endometrialIds = set(impactSigs[impactSigs['cancer_type'] == 'Endometrial Cancer']['Tumor_Sample_Barcode'])
endometrialHypermutatorIds = set(impactSigs[(impactSigs['cancer_type'] == 'Endometrial Cancer') & (impactSigs['Nmut_Mb'] > N_MUT_MB_THRESHOLD)]['Tumor_Sample_Barcode'])
endometrialMuts = mafWithHotspotInfo[mafWithHotspotInfo['Tumor_Sample_Barcode'].isin(endometrialIds)]
endometrialSummaryDf = summarize_n_double_mutations(endometrialMuts, sigsNmutDict)
endometrialSummaryDf['isHypermutator'] = endometrialSummaryDf['Tumor_Sample_Barcode'].apply(lambda x: True if x in endometrialHypermutatorIds else False)
#endometrialHyperMuts = mafWithHotspotInfo[mafWithHotspotInfo['Tumor_Sample_Barcode'].isin(endometrialHypermutatorIds)]
endometrialHypermutationSummaryDf = summarize_n_double_mutations(endometrialHyperMuts, sigsNmutDict)
endometrialHypermutationSummaryDf[endometrialHypermutationSummaryDf['nGenesWithDoubletHotspotMutations'] > 0].shape[0]*1.0 /endometrialHypermutationSummaryDf.shape[0]
endometrialHypermutationSummaryDf['cancerType'] = 'Endometrial Cancer'

#########################################################
#########APOBEC ANALYSES

bladderHypermutatorIds = set(impactSigs[(impactSigs['cancer_type'] == 'Bladder Cancer') & (impactSigs['Nmut_Mb'] > N_MUT_MB_THRESHOLD)]['Tumor_Sample_Barcode'])
bladderHypermutatorMuts = mafWithHotspotInfo[mafWithHotspotInfo['Tumor_Sample_Barcode'].isin(bladderHypermutatorIds)]
bladderHypermutationSummaryDf = summarize_n_double_mutations(bladderHypermutatorMuts, sigsNmutDict)
bladderHypermutationSummaryDf['cancerType'] = 'Bladder Cancer'

###############################################
############COLORECTAL ANALYSES
colorectalHypermutatorIds = set(impactSigs[(impactSigs['cancer_type'] == 'Colorectal Cancer') & (impactSigs['Nmut_Mb'] > N_MUT_MB_THRESHOLD)]['Tumor_Sample_Barcode'])
colorectalHyperMuts = mafWithHotspotInfo[mafWithHotspotInfo['Tumor_Sample_Barcode'].isin(colorectalHypermutatorIds)]
colorectalHypermutationSummaryDf = summarize_n_double_mutations(colorectalHyperMuts, sigsNmutDict)
colorectalHypermutationSummaryDf['cancerType'] = 'Colorectal Cancer'

##########################################
#############TMZ Gliomas 
gliomaHypermutatorTMZIds = set(impactSigs[(impactSigs['cancer_type'] == 'Glioma') & (impactSigs['Nmut_Mb'] > N_MUT_MB_THRESHOLD)]['Tumor_Sample_Barcode'])
gliomaHypermutatorMuts = mafWithHotspotInfo[mafWithHotspotInfo['Tumor_Sample_Barcode'].isin(gliomaHypermutatorTMZIds)]
gliomaHypermutationSummaryDf = summarize_n_double_mutations(gliomaHypermutatorMuts, sigsNmutDict)
gliomaHypermutationSummaryDf['cancerType'] = 'Glioma'

#########################################
############MELANOMAS
melanomaHypermutatorIds = set(impactSigs[(impactSigs['cancer_type'] == 'Melanoma') & (impactSigs['Nmut_Mb'] > N_MUT_MB_THRESHOLD)]['Tumor_Sample_Barcode'])
melanomaHypermutatorMuts = mafWithHotspotInfo[mafWithHotspotInfo['Tumor_Sample_Barcode'].isin(melanomaHypermutatorIds)]
melanomaHypermutationSummaryDf = summarize_n_double_mutations(melanomaHypermutatorMuts, sigsNmutDict)
melanomaHypermutationSummaryDf['cancerType'] = 'Melanoma'

#################################3
###############BREAST CANCERS
breastHypermutatorIds = set(impactSigs[(impactSigs['cancer_type'] == 'Breast Cancer') & (impactSigs['Nmut_Mb'] > N_MUT_MB_THRESHOLD)]['Tumor_Sample_Barcode'])
breastHypermutatorMuts = mafWithHotspotInfo[mafWithHotspotInfo['Tumor_Sample_Barcode'].isin(breastHypermutatorIds)]
breastHypermutationSummaryDf = summarize_n_double_mutations(breastHypermutatorMuts, sigsNmutDict)
breastHypermutationSummaryDf['cancerType'] = 'Breast Cancer'"""

mutsInCancerTypesOfInterest = mafWithHotspotInfo[(mafWithHotspotInfo['cancer_type'] == 'Breast Cancer') |
        (mafWithHotspotInfo['cancer_type'] == 'Endometrial Cancer') |
        (mafWithHotspotInfo['cancer_type'] == 'Bladder Cancer') |
        (mafWithHotspotInfo['cancer_type'] == 'Colorectal Cancer') |
        (mafWithHotspotInfo['cancer_type'] == 'Glioma') |
        (mafWithHotspotInfo['cancer_type'] == 'Melanoma') ]
    
mutsInCancerTypesOfInterest['isHypermutator'] = mutsInCancerTypesOfInterest['Tumor_Sample_Barcode'].apply(lambda x: True if x in sigsNmutDict and sigsNmutDict[x] >  N_MUT_MB_THRESHOLD else False)
dfSummary = summarize_n_double_mutations(mutsInCancerTypesOfInterest, sigsNmutDict)
dfCombined = dfSummary
dfCombined['pid'] = dfCombined['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
dfCombined['cancerType'] = dfCombined['pid'].apply(lambda x: cDict[x] if x in cDict else None)

#dfCombined = pd.concat([endometrialHypermutationSummaryDf, bladderHypermutationSummaryDf, colorectalHypermutationSummaryDf, gliomaHypermutationSummaryDf, melanomaHypermutationSummaryDf, breastHypermutationSummaryDf])
dfCombined['doubleHotspotPresent'] = dfCombined['nGenesWithDoubletHotspotMutations'].apply(lambda x: 1 if x > 0 else 0)
dfCombined['signatureType'] = dfCombined['Tumor_Sample_Barcode'].apply(lambda x:
    'POLE' if x in poleHypermutatorIds
    else 'MMR' if x in mmrHypermutatorIds
    else 'APOBEC' if x in apobecHypermutatorIds
    else 'TMZ' if x in tmzHypermutatorIds
    else 'UV' if x in uvHypermutatorIds
    else 'Other')
dfCombined['hotspotsNames'] = dfCombined['genesWithDoubletHotspots'].apply(lambda x: ';'.join(x))
dfCombined['hotspotMutated'] = dfCombined['hotspotsNames'].apply(lambda x:
    'PIK3CA' if x == 'PIK3CA'
    else 'TP53' if x == 'TP53'
    else 'PTEN' if x == 'PTEN'
    else 'ERBB2' if x == 'ERB2'
    else 'multiple(inc:PIK3CA)' if (';' in x and 'PIK3CA' in x)
    else 'multiple' if ';' in x
    else 'other' if len(x) > 0
    else 'no doublets'
    )
    
dfCombined['cancerSubtype'] = dfCombined['cancerType']
dfCombined['cancerSubtype'] = dfCombined.apply(lambda row:
    'EndometrialPOLE' if row['cancerType'] == 'Endometrial Cancer' and row['signatureType'] == 'POLE'
    else 'EndometrialMMR' if row['cancerType'] == 'Endometrial Cancer'
    else 'ColorectalPOLE 'if row['cancerType'] == 'Colorectal Cancer' and row['signatureType'] == 'POLE'
    else 'ColorectalMMR 'if row['cancerType'] == 'Colorectal Cancer'
    else row['cancerType']
,axis=1)
    
dfCombined['cancerAndMutationSubtype'] = dfCombined.apply(lambda row: row['cancerType'] + '_Hypermutator' if row['isHypermutator'] else row['cancerType'] + '_nonHypermutator', axis=1)
    
dfCombined.to_csv('~/Desktop/dataForLocalPlotting/hypermutatorDoubletAnalysis.tsv', sep='\t', index=False)


