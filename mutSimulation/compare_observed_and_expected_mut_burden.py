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
#import maf_analysis_utils
import mutation_modeling_util

#we look at the size of the impact panel to determine what fraction of the panel are tumor suppresors
#TODO create an indel oncogenicity for each panel
def derive_indel_oncogenicity():
    fractionOfIndelsThatAreFrameShift = .9 #calculated with craig (approximately 90% of msi-indels in oncogenes are frameshifts)
    repeatRegionInfo = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/IMPACT_repeat_stats.txt')
    tumorSuppressors = set(['ERRFI1', 'ASXL2', 'PMAIP1', 'ACTG1', 'SUFU', 'FBXO11', 'MEN1', 'FAM58A', 'B2M', 'RB1', 'DUSP22', 'SESN1', 'GPS2', 'RAD51D', 'SMG1', 'CDC73', 'MAP3K1', 'SMARCB1', 'INPP4B', 'PARK2', 'SMAD4', 'CBFB', 'CDH1', 'PPP6C', 'SETDB1', 'SETDB2', 'NF2', 'CDKN2B', 'CDKN2C', 'CDKN2A', 'DDX3X', 'PIK3R1', 'BARD1', 'PDS5B', 'KLF4', 'SPRED1', 'VHL', 'SMAD2', 'PMS1', 'PMS2', 'SETD2', 'GATA3', 'TBL1XR1', 'MUTYH', 'SOCS1', 'FAM175A', 'ROBO1', 'ARID1B', 'ARID1A', 'TCF7L2', 'STK11', 'FOXA1', 'PTEN', 'FAT1', 'FAS', 'CYLD', 'MAX', 'SH2D1A', 'APC', 'NTHL1', 'CTCF', 'KDM5C', 'KMT2C', 'ZFHX3', 'FOXP1', 'PIGA', 'CDKN1B', 'CDKN1A', 'FUBP1', 'MSH2', 'ID3', 'TNFRSF14', 'TRAF3', 'EP400', 'BRIP1', 'ARID4A', 'ARID4B', 'XRCC2', 'DAXX', 'SDHAF2', 'ASXL1', 'AMER1', 'RASA1', 'EGR1', 'MST1', 'SOX17', 'RUNX1', 'PIK3R3', 'NCOR1', 'NF1', 'JAK1', 'PTPRD', 'CHEK2', 'CHEK1', 'SMC1A', 'TMEM127', 'STAG1', 'RAD51', 'TCF3', 'STAG2', 'ARID2', 'RAD50', 'RNF43', 'PARP1', 'BLM', 'CUX1', 'RECQL', 'RAD21', 'PTPN2', 'PTPN1', 'SLX4', 'INHA', 'PAX5', 'IRF1', 'TP53', 'HLA-A', 'IRF8', 'CBL', 'TOP1', 'SHQ1', 'PRDM1', 'NSD1', 'ATXN2', 'CREBBP', 'HDAC4', 'SESN2', 'PPP2R1A', 'EPHA7', 'ATM', 'EPHA3', 'POT1', 'SMAD3', 'MOB3B', 'TBX3', 'POLE', 'ATR', 'FANCD2', 'FH', 'BCORL1', 'SOX9', 'IKZF3', 'TSC1', 'TP63', 'MRE11A', 'SDHC', 'BTG1', 'POLD1', 'CIITA', 'SMC3', 'SAMHD1', 'RTEL1', 'ECT2L', 'PIK3R2', 'CRBN', 'FANCC', 'NBN', 'FANCA', 'HLA-B', 'RECQL4', 'DUSP4', 'ERCC2', 'FBXW7', 'TGFBR2', 'TGFBR1', 'MSH3', 'RBM15', 'TET1', 'TET3', 'SESN3', 'MGA', 'LTB', 'FOXL2', 'SH2B3', 'BCOR', 'HIST1H1D', 'ATRX', 'EP300', 'RAD51C', 'RAD51B', 'HIST1H1B', 'TNFAIP3', 'DICER1', 'ARID5B', 'LATS2', 'FOXO1', 'KEAP1', 'EZH2', 'SP140', 'NKX3-1', 'PBRM1', 'PALB2', 'CIC', 'BRCA1', 'DTX1', 'FLCN', 'SPEN', 'CD58', 'ERCC3', 'ERCC4', 'MSH6', 'BCL11B', 'BMPR1A', 'ERF', 'BRCA2', 'NOTCH2', 'EED', 'MITF', 'ELF3', 'SMARCA4', 'BBC3', 'ANKRD11', 'CEBPA', 'BCL2L11', 'AXIN2', 'AXIN1', 'CDK12', 'ESCO2', 'MLH1', 'SDHB', 'MED12', 'HNF1A', 'RYBP', 'ATP6V1B2', 'DNMT3B', 'KMT2B', 'KMT2A', 'DNMT3A', 'NFKBIA', 'TRAF5', 'KMT2D', 'SPOP', 'RBM10', 'P2RY8', 'TP53BP1', 'TSC2', 'KDM6A', 'EPCAM', 'PHOX2B', 'NPM1', 'BCL10', 'LATS1', 'HOXB13', 'ARID3A', 'PTPRT', 'PTPRS', 'INPPL1', 'NOTCH4', 'TET2', 'NOTCH1', 'CASP8', 'NOTCH3', 'GRIN2A', 'MAP2K4', 'WT1', 'BACH2', 'SDHA', 'BAP1', 'PTCH1', 'SDHD'])
    return fractionOfIndelsThatAreFrameShift*np.nansum(repeatRegionInfo[repeatRegionInfo['Hugo_Symbol'].isin(tumorSuppressors)]['bp'])/np.nansum(repeatRegionInfo['bp'])

#little function that loads indel data
def get_indel_dict():
    returnDict = dict()
    sigs = ['Signature.' + str(i) for i in range(1,31)]
    indelFracDf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/indelRateBySignature.tsv')
    indelFracDict = dict(zip(indelFracDf['Signature'], indelFracDf['Indel_Fraction']))
    for sig in sigs:
        if sig in indelFracDict:
            returnDict[sig] = indelFracDict[sig]
        else:
            returnDict[sig] = indelFracDict['All_IMPACT_Cases']
    return returnDict

def get_per_case_mut_info(nmutDfPath = '/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/nmutInfo_impact_filtered.tsv'):
    df = pd.read_table(nmutDfPath)
    return dict(zip(df['Tumor_Sample_Barcode'], df['Nmut']))

def get_per_case_oncogenic_mut_info(muts):
    oncogenicMuts = muts[muts['oncogenic'].notnull()]
    nMutOncDict = dict(oncogenicMuts['Tumor_Sample_Barcode'].value_counts())
    return nMutOncDict

#TODO flesh this out
def get_per_case_unique_gene_oncogenic_mut_info(muts):
    oncogenicMuts = muts[muts['oncogenic'].notnull()]
    oncogenicMuts['caseGeneAltered'] = oncogenicMuts.apply(lambda row: row['Tumor_Sample_Barcode'] + '_' + row['Hugo_Symbol'], axis=1)
    oncogenicMuts = oncogenicMuts.drop_duplicates(subset=['caseGeneAltered'])
    return dict(oncogenicMuts['Tumor_Sample_Barcode'].value_counts())

def get_per_gene_oncogenic_mut_info(muts):
    oncogenicMuts = muts[muts['oncogenic'].notnull()]
    nMutOncGeneDict = dict(oncogenicMuts['Hugo_Symbol'].value_counts())
    return nMutOncGeneDict

def get_per_case_indel_fracs(muts):
    d = {}
    indelClassificationNames = set(['Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins'])
    nMutDict = dict(muts['Tumor_Sample_Barcode'].value_counts())
    nIndelDict = dict(muts[muts['Variant_Classification'].isin(indelClassificationNames)]['Tumor_Sample_Barcode'].value_counts())
    cases = set(muts['Tumor_Sample_Barcode'])
    for case in cases:
        nMut = 0
        nIndel = 0
        if case in nMutDict:
            nMut = nMutDict[case]
        if case in nIndelDict:
            nIndel = nIndelDict[case]
        if nMut == 0:
            d[case] = 0
        else:
            d[case] = 1.0*nIndel/nMut
    return d

def fix_mll_genes(maf):
    maf['Hugo_Symbol'] = maf['Hugo_Symbol'].apply(lambda x:
        'KMT2A' if x == 'MLL'
        else 'KMT2B' if x == 'MLL2'
        else 'KMT2C' if x == 'MLL3'
        else x)   
    return maf

#EXECUTABLE FROM THE COMMAND LINE CAUSE IT TAKES FOREVER
#REAL data
print 'running '

casesFile = sys.argv[1]
mode = sys.argv[2]
outfile = sys.argv[3]

f = open(casesFile)
cases = f.readlines()
cases = set([i.strip('\n') for i in cases])

print 'readfing in signatures'
impactSigs = pd.read_table(pathPrefix + '/ifs/res/taylorlab/impact_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
print 'reading in mutations maf'
impactMutsMaf = analysis_utils.load_in_df_with_progress(pathPrefix + '/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/all_impact_mutations_annotated_cohort.maf', 275171)

mutsToConsider = impactMutsMaf[impactMutsMaf['Tumor_Sample_Barcode'].isin(cases)]
mutsToConsider = fix_mll_genes(mutsToConsider)
sigsToConsider = impactSigs[impactSigs['Tumor_Sample_Barcode'].isin(cases)]

print 'summarizing oncogenic observed and indels'
#NOW START GETTING DATA
indelProbDict = get_indel_dict()
#indelOncogenicity = derive_indel_oncogenicity()
perCaseIndelFracDict = get_per_case_indel_fracs(mutsToConsider)
nMutDict = get_per_case_mut_info()
perCaseOncogenicInfo = get_per_case_oncogenic_mut_info(mutsToConsider)
perCaseUniqueOncogenicInfo = get_per_case_unique_gene_oncogenic_mut_info(mutsToConsider)

print 'reading in sim data'
simOncogenicitySummary = pd.read_table('/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/simulatedMutationSummary.tsv')

if mode == 'perCase':

    oncogenicSusceptibilityDict = mutation_modeling_util.calculate_quadnuc_based_oncogenic_susceptibility_dict(simOncogenicitySummary)
    print 'comparing observed and expected'

    #WE add the MSI class to the impact signatures
    impactSigs = analysis_utils.map_cases_to_msi_sensor_class(impactSigs, msiSensorInfo=pathPrefix + '/ifs/work/taylorlab/friedman/mskImpactAsOfMarch2019/dmp/mskimpact/data_clinical_sample.txt')

    df = mutation_modeling_util.get_observed_expected_values_for_cohort(cases, impactSigs, oncogenicSusceptibilityDict, nMutDict, perCaseOncogenicInfo, perCaseUniqueOncogenicInfo, perCaseIndelFracDict)

    filePath = '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/observedVsExpectedData.tsv'
    print 'writing file to ', filePath
    df.to_csv(filePath, index=False, sep='\t')

elif mode == 'perGene':
    oncogenicMutsObservedPerGeneDict = get_per_gene_oncogenic_mut_info(mutsToConsider)
    oncogenicSusceptibilityDictByGene = mutation_modeling_util.calculate_quadnuc_and_gene_based_oncogenic_susceptibility_dict(simOncogenicitySummary)

    print 'comparing observed and expected'
    df = mutation_modeling_util.get_observed_expected_values_for_genes(cases, #the tumor sample barcodes of cases we ought to consider
        nMutDict,
        impactSigs, #the decomposed signatures of the cases we ought to consider
        oncogenicSusceptibilityDictByGene, #a dictionary mapping quadnuc muts in each gene to oncogenic probabilities
        oncogenicMutsObservedPerGeneDict #a dictionary that tells us how many oncogenic mutations we observe in a gene in a cohort
    )

    print 'writing file to ', outfile
    df.to_csv(outfile, index=False, sep='\t')



else: print 'improper mode specified'


#mutSimulationIds.txt
#usage python compare_observed_and_expected_mut_burden.py ~/friedman/myAdjustedDataFiles/mutSimulationIdsSmall.txt perCase
#python compare_observed_and_expected_mut_burden.py /ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/idFiles/endometrialHypermutantIds.txt perGene


