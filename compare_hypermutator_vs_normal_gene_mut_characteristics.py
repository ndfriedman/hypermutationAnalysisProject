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

def summarize_per_case_uniq_muts_info(mafX, recurrentGenes):
    oncMutUniq = []
    oncMutRecurrentUniq = []
    for case in set(mafX['Tumor_Sample_Barcode']):
        caseMaf = mafX[mafX['Tumor_Sample_Barcode'] == case]
        caseMafRecurrent = caseMaf[caseMaf['Hugo_Symbol'].isin(recurrentGenes)]
        caseMafNotRecurrent = caseMaf[~caseMaf['Hugo_Symbol'].isin(recurrentGenes)]
        nGenesOncogenicallyMutated = len(set(caseMafNotRecurrent['Hugo_Symbol']))
        nRecurrentGenesOncogenicallyMutated = len(set(caseMafRecurrent['Hugo_Symbol']))
        oncMutUniq.append(nGenesOncogenicallyMutated)
        oncMutRecurrentUniq.append(nRecurrentGenesOncogenicallyMutated)
    return np.nanmean(oncMutUniq), np.nanmean(oncMutRecurrentUniq)

#Note the maf passed in here should be oncogenic muts only
def count_case_and_cohort_gene_mut_info(mafY, recurrentGenes):
    nCases = len(set(mafY['Tumor_Sample_Barcode']))
    nRecurrentMutsPerCase = 1.0*mafY[mafY['Hugo_Symbol'].isin(recurrentGenes)].shape[0]/nCases
    nNotRecurrentMutsPerCase =  1.0*mafY[~mafY['Hugo_Symbol'].isin(recurrentGenes)].shape[0]/nCases
    oncMutUniq, oncMutRecurrentUniq = summarize_per_case_uniq_muts_info(mafY, recurrentGenes)
    return nRecurrentMutsPerCase, nNotRecurrentMutsPerCase, oncMutUniq, oncMutRecurrentUniq
    
def summarize_data_across_cancer_types(maf, cancerTypes, sigs, notHypermutationThresh = 20, hypermutationThresh = 50):
    listOfDicts = []
    oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
    for cancerType in cancerTypes:
        cancerTypeMaf = maf[maf['cancer_type'] == cancerType]
        cancerTypeHyperIds = set(sigs[(sigs['cancer_type'] == cancerType) & (sigs['Nmut_Mb'] > hypermutationThresh)]['Tumor_Sample_Barcode'])
        cancerTypeNormalIds = set(sigs[(sigs['cancer_type'] == cancerType) & (sigs['Nmut_Mb'] < notHypermutationThresh)]['Tumor_Sample_Barcode'])
        recurrentTumorSupressors, recurrentOncogenes = maf_analysis_utils.enumerate_recurrently_mutated_tumor_supressors_and_oncogenes(cancerTypeMaf[cancerTypeMaf['Tumor_Sample_Barcode'].isin(cancerTypeNormalIds)], thresh=.05)
        oncogenicHyperMaf = cancerTypeMaf[(cancerTypeMaf['oncogenic'].isin(oncogenicMutColNames)) & (cancerTypeMaf['Tumor_Sample_Barcode'].isin(cancerTypeHyperIds))]
        normalMaf = cancerTypeMaf[(cancerTypeMaf['oncogenic'].isin(oncogenicMutColNames)) & (cancerTypeMaf['Tumor_Sample_Barcode'].isin(cancerTypeNormalIds))]
        
        #Get information and put into a dictionary to return
        #HYPERMUTATED
        print 'Analyzing', cancerType
        nRecurrentMutsPerCase, nNotRecurrentMutsPerCase, oncMutUniq, oncMutRecurrentUniq = count_case_and_cohort_gene_mut_info(oncogenicHyperMaf, recurrentTumorSupressors | recurrentOncogenes)
        listOfDicts.append({'class': cancerType + 'Related Genes', 'Nmut_case': nRecurrentMutsPerCase, 'NUniqMut_case': oncMutRecurrentUniq, 'Hypermutated': True})
        listOfDicts.append({'class': cancerType + ' (other)', 'Nmut_case': nNotRecurrentMutsPerCase, 'NUniqMut_case': oncMutUniq, 'Hypermutated': True})
        
        #NOT HYPERMUTATED
        nRecurrentMutsPerCase, nNotRecurrentMutsPerCase, oncMutUniq, oncMutRecurrentUniq = count_case_and_cohort_gene_mut_info(normalMaf, recurrentTumorSupressors | recurrentOncogenes)
        listOfDicts.append({'class': cancerType + ' Related Genes', 'Nmut_case': nRecurrentMutsPerCase, 'NUniqMut_case': oncMutRecurrentUniq, 'Hypermutated': False})
        listOfDicts.append({'class': cancerType + ' (other)', 'Nmut_case': nNotRecurrentMutsPerCase, 'NUniqMut_case': oncMutUniq, 'Hypermutated': False})
        
    return pd.DataFrame(listOfDicts)

#TODO ADD TUMOR SUPPRESOR FUNCTIONALITY

def enumerate_related_and_unrelated_oncogenic_mut_info(fullMaf, nmutMbInfo, cancerTypeDict, recurrentGenesDict):
    listOfDicts = []
    oncogenicMaf = fullMaf[fullMaf['oncogenic'].notnull()]
    cntr = 0
    nCases = len(set(oncogenicMaf['Tumor_Sample_Barcode']))
    for case in set(oncogenicMaf['Tumor_Sample_Barcode']):
        cntr += 1
        if cntr%100 == 0: print cntr, nCases
        if case in cancerTypeDict and case in nmutMbInfo:
            cType = cancerTypeDict[case]
            if cType in recurrentGenesDict.keys():
                caseMaf = fullMaf[fullMaf['Tumor_Sample_Barcode'] == case]
                
                #DATA before consideration of tumor suppressor vs oncogene
                caseMafCTypeRelated = caseMaf[caseMaf['Hugo_Symbol'].isin(recurrentGenesDict[cType])]
                caseMafCTypeNotRelated = caseMaf[~caseMaf['Hugo_Symbol'].isin(recurrentGenesDict[cType])]
                caseMafCTypeRelatedUniq = caseMafCTypeRelated.drop_duplicates(subset=['Hugo_Symbol'])
                caseMafCTypeNotRelatedUniq = caseMafCTypeNotRelated.drop_duplicates(subset=['Hugo_Symbol'])

                #PLOT tumor suppressor/oncogene info                
                tumorSupressors = set(['ERRFI1', 'ASXL2', 'PMAIP1', 'ACTG1', 'SUFU', 'FBXO11', 'MEN1', 'FAM58A', 'B2M', 'RB1', 'DUSP22', 'SESN1', 'GPS2', 'RAD51D', 'SMG1', 'CDC73', 'MAP3K1', 'SMARCB1', 'INPP4B', 'PARK2', 'SMAD4', 'CBFB', 'CDH1', 'PPP6C', 'SETDB1', 'SETDB2', 'NF2', 'CDKN2B', 'CDKN2C', 'CDKN2A', 'DDX3X', 'PIK3R1', 'BARD1', 'PDS5B', 'KLF4', 'SPRED1', 'VHL', 'SMAD2', 'PMS1', 'PMS2', 'SETD2', 'GATA3', 'TBL1XR1', 'MUTYH', 'SOCS1', 'FAM175A', 'ROBO1', 'ARID1B', 'ARID1A', 'TCF7L2', 'STK11', 'FOXA1', 'PTEN', 'FAT1', 'FAS', 'CYLD', 'MAX', 'SH2D1A', 'APC', 'NTHL1', 'CTCF', 'KDM5C', 'KMT2C', 'ZFHX3', 'FOXP1', 'PIGA', 'CDKN1B', 'CDKN1A', 'FUBP1', 'MSH2', 'ID3', 'TNFRSF14', 'TRAF3', 'EP400', 'BRIP1', 'ARID4A', 'ARID4B', 'XRCC2', 'DAXX', 'SDHAF2', 'ASXL1', 'AMER1', 'RASA1', 'EGR1', 'MST1', 'SOX17', 'RUNX1', 'PIK3R3', 'NCOR1', 'NF1', 'JAK1', 'PTPRD', 'CHEK2', 'CHEK1', 'SMC1A', 'TMEM127', 'STAG1', 'RAD51', 'TCF3', 'STAG2', 'ARID2', 'RAD50', 'RNF43', 'PARP1', 'BLM', 'CUX1', 'RECQL', 'RAD21', 'PTPN2', 'PTPN1', 'SLX4', 'INHA', 'PAX5', 'IRF1', 'TP53', 'HLA-A', 'IRF8', 'CBL', 'TOP1', 'SHQ1', 'PRDM1', 'NSD1', 'ATXN2', 'CREBBP', 'HDAC4', 'SESN2', 'PPP2R1A', 'EPHA7', 'ATM', 'EPHA3', 'POT1', 'SMAD3', 'MOB3B', 'TBX3', 'POLE', 'ATR', 'FANCD2', 'FH', 'BCORL1', 'SOX9', 'IKZF3', 'TSC1', 'TP63', 'MRE11A', 'SDHC', 'BTG1', 'POLD1', 'CIITA', 'SMC3', 'SAMHD1', 'RTEL1', 'ECT2L', 'PIK3R2', 'CRBN', 'FANCC', 'NBN', 'FANCA', 'HLA-B', 'RECQL4', 'DUSP4', 'ERCC2', 'FBXW7', 'TGFBR2', 'TGFBR1', 'MSH3', 'RBM15', 'TET1', 'TET3', 'SESN3', 'MGA', 'LTB', 'FOXL2', 'SH2B3', 'BCOR', 'HIST1H1D', 'ATRX', 'EP300', 'RAD51C', 'RAD51B', 'HIST1H1B', 'TNFAIP3', 'DICER1', 'ARID5B', 'LATS2', 'FOXO1', 'KEAP1', 'EZH2', 'SP140', 'NKX3-1', 'PBRM1', 'PALB2', 'CIC', 'BRCA1', 'DTX1', 'FLCN', 'SPEN', 'CD58', 'ERCC3', 'ERCC4', 'MSH6', 'BCL11B', 'BMPR1A', 'ERF', 'BRCA2', 'NOTCH2', 'EED', 'MITF', 'ELF3', 'SMARCA4', 'BBC3', 'ANKRD11', 'CEBPA', 'BCL2L11', 'AXIN2', 'AXIN1', 'CDK12', 'ESCO2', 'MLH1', 'SDHB', 'MED12', 'HNF1A', 'RYBP', 'ATP6V1B2', 'DNMT3B', 'KMT2B', 'KMT2A', 'DNMT3A', 'NFKBIA', 'TRAF5', 'KMT2D', 'SPOP', 'RBM10', 'P2RY8', 'TP53BP1', 'TSC2', 'KDM6A', 'EPCAM', 'PHOX2B', 'NPM1', 'BCL10', 'LATS1', 'HOXB13', 'ARID3A', 'PTPRT', 'PTPRS', 'INPPL1', 'NOTCH4', 'TET2', 'NOTCH1', 'CASP8', 'NOTCH3', 'GRIN2A', 'MAP2K4', 'WT1', 'BACH2', 'SDHA', 'BAP1', 'PTCH1', 'SDHD'])
                caseMafCTypeRelatedTumorSuppressors = caseMafCTypeRelated[caseMafCTypeRelated['Hugo_Symbol'].isin(tumorSupressors)]
                caseMafCTypeRelatedOncogenes = caseMafCTypeRelated[~caseMafCTypeRelated['Hugo_Symbol'].isin(tumorSupressors)]
                caseMafCTypeNotRelatedTumorSuppressors = caseMafCTypeNotRelated[caseMafCTypeNotRelated['Hugo_Symbol'].isin(tumorSupressors)]
                caseMafCTypeNotRelatedOncogenes = caseMafCTypeNotRelated[~caseMafCTypeNotRelated['Hugo_Symbol'].isin(tumorSupressors)]
                
                listOfDicts.append(
                        {'Nmut_Related': caseMafCTypeRelated.shape[0], 
                         'Nmut_Non-Related': caseMafCTypeNotRelated.shape[0],
                         'Nmut_Uniq_Related': caseMafCTypeRelatedUniq.shape[0], 
                         'Nmut_Uniq_Non-Related': caseMafCTypeNotRelatedUniq.shape[0],
                         
                         'Nmut_Related_Oncogene': caseMafCTypeRelatedOncogenes.shape[0], 
                         'Nmut_Non-Related_Oncogene': caseMafCTypeNotRelatedOncogenes.shape[0],
                         'Nmut_Related_Tumor_Supressor': caseMafCTypeRelatedTumorSuppressors.shape[0], 
                         'Nmut_Non-Related_Tumor_Supressor': caseMafCTypeNotRelatedTumorSuppressors.shape[0],
                         
                         'Nmut_Mb_Case': nmutMbInfo[case],
                         'Tumor_Sample_Barcode': case
                })
    return pd.DataFrame(listOfDicts)

filteredMafDf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/tempScriptFiles/filteredMafAnnotated.maf')
filteredMafDf['pid'] = filteredMafDf['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
filteredMafDf = maf_analysis_utils.fix_mll_genes(filteredMafDf)

#Get cancer type and nmut mb info from signatures
sigsInfo = pd.read_table(pathPrefix + '/ifs/res/taylorlab/impact_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
sigsInfo['pid'] = sigsInfo['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/mskImpactAsOfMarch2019/dmp/mskimpact/data_clinical_sample.txt')
sigsInfo['cancer_type'] = sigsInfo['pid'].apply(lambda x: cDict[x] if x in cDict else None)

nmutMbDict = dict(zip(sigsInfo['Tumor_Sample_Barcode'], sigsInfo['Nmut_Mb']))
cTypeDict = dict(zip(sigsInfo['Tumor_Sample_Barcode'], sigsInfo['cancer_type']))

cTypes = set([re.sub('_', ' ', x.strip('.tsv')) for x in os.listdir(pathPrefix + '/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/hypermutationStatusIds')])
cTypes.remove('Cancer of Unknown Primary')
mutatedGenes = maf_analysis_utils.create_dictionary_mapping_genes_to_cancer_types_with_implication(filteredMafDf, cancerTypes = cTypes)

df = enumerate_related_and_unrelated_oncogenic_mut_info(filteredMafDf, nmutMbDict, cTypeDict, mutatedGenes)
df = analysis_utils.map_cases_to_msi_sensor_class(df, msiSensorInfo = pathPrefix + '/ifs/work/taylorlab/friedman/mskImpactAsOfMarch2019/dmp/mskimpact/data_clinical_sample.txt')

df.to_csv('~/Desktop/WORK/dataForLocalPlotting/oncogenicMutCount.tsv', sep='\t', index=False)


d = maf_analysis_utils.calculate_nmut_mb_info_from_filtered_maf(filteredMafDf)
d['pid'] = d['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
d['cancer_type'] = d['pid'].apply(lambda x: cDict[x] if x in cDict else None)

d.to_csv('~/Desktop/mnt/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/tmbInfo.tsv', index=False, sep='\t')









###############
    #################3
    #################
#OLD VERSION TO DELETE???

#mafWithInfo = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/annotatedOncoPlusHotspotMafAllImpact_trinuc')
mafWithInfo = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/tempScriptFiles/attempt2Maf.maf')
impactSigs = pd.read_table(pathPrefix + '/ifs/res/taylorlab/impact_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
impactSigs['pid'] = impactSigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
impactSigs['cancer_type'] = impactSigs['pid'].apply(lambda x: cDict[x] if x in cDict else None)
mafWithInfo['pid'] = mafWithInfo['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
mafWithInfo['cancer_type'] = mafWithInfo['pid'].apply(lambda x: cDict[x] if x in cDict else None)

#adjust column names to make the 
renameDict = {key:value for (key,value) in [('mean_' + str(i), 'Signature.' + str(i)) for i in range(1,31)]}
impactSigs = impactSigs.rename(columns=renameDict)
impactSigs = mutationSigUtils.merge_signature_columns(impactSigs, mode='Stratton', drop=False, smokingMerge=True, confidence=False, mean=True, prefix='Signature.')
impactSigs = mutationSigUtils.merge_signature_columns(impactSigs, mode='Stratton', drop=True, smokingMerge=False, confidence=True, mean=True, prefix='mean_')

df = summarize_data_across_cancer_types(mafWithInfo, ['Endometrial Cancer', 'Colorectal Cancer', 'Glioma'], impactSigs)
df['orderingVal'] = df['class'].apply(lambda x: 1 if 'Endometrial' in x else 2 if 'Colorectal' in x else 3)
df['orderingVal'] = df.apply(lambda row: row['orderingVal'] + 0.1 if 'Not' in row['class'] else row['orderingVal'], axis=1)
df.to_csv('~/Desktop/WORK/dataForLocalPlotting/hyperVsNormalMutCharComp.tsv', index=False, sep='\t')


mafWithInfo = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/annotatedOncoPlusHotspotMafAllImpact_trinuc')


##########3TEMPORARY ANALYSIS





