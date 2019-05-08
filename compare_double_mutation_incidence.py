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

import scipy.stats

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h


#returns proto dataframe information with averages of each column
def process_mut_summary_maf_for_averages(summaryMaf, genes, suffix='_doubleOncogenic', labelMode='hypermutator'):
    
    listOfDicts = [] #the list of dictionaries we will return from this function
    sortedArr, allValues, avg = summarize_top_double_mutated_genes(summaryMaf, suffix)
    m, lowerBound, upperBound = mean_confidence_interval(allValues, confidence=0.95)
    listOfDicts.append({'label': labelMode + '_all_genes_average', 'value': avg, 'lowerConf':lowerBound, 'upperConf':upperBound})
    for gene in genes:
        col = gene + suffix
        avg = np.nanmean(summaryMaf[col])
        m, lowerBound, upperBound = mean_confidence_interval(list(summaryMaf[col]), confidence=0.95)
        listOfDicts.append({'label': labelMode + '_' + gene, 'value': avg, 'lowerConf':lowerBound, 'upperConf':upperBound})
    return listOfDicts

#a utility function that summarizes the genes that are most recurrenctly double mutated across a cohort
#returns a sorted list of the genes, all 1s/0s for what happens for each level, the average
def summarize_top_double_mutated_genes(sMaf, suffix='_doubleOncogenic'):
    tupleArr = []
    valueArr = []
    arr = []
    for column in sMaf.columns.values:
        if suffix in column:
            arr.append(column)
    for val in arr: 
        tupleArr.append((val, np.nanmean(sMaf[val])))
        valueArr = valueArr + list(sMaf[val])
    return sorted(tupleArr, key= lambda x: x[1], reverse=True), valueArr, np.nanmean([i[1] for i in tupleArr])

def assign_ordering_val(row, orderingDict):
    geneName = row['label'].split('_')[len(row['label'].split('_')) - 1] #get the substring after the last underscore (the gene name or the word "average")
    isHypermutator = True
    if row['label'][0] == 'N':
        isHypermutator = False  #the string starts with NOT-hypermutator
    if geneName == 'average':
        if isHypermutator:
            return 1000 #big value to make sure these stay at the end
        else:
            return 1001
    else:
        if isHypermutator:
            return orderingDict[geneName]
        else:
            return orderingDict[geneName] + .1

#big function that takes a cohort maf and generates statistics on doublet prevalence of nDoublesToPlot doublets     
def create_hypermutation_prevalence_summary(mutSummaryMaf, nDoublesToPlot = 10, doubleSuffix='_doubleOncogenic'):
    mutSummaryMafHypermutators = mutSummaryMaf[mutSummaryMaf['isHypermutator'] == True]
    mutSummaryMafNonHypermutators = mutSummaryMaf[mutSummaryMaf['isHypermutator'] == False]
    
    topHypermutatorGenesArr, v1, v2 = summarize_top_double_mutated_genes(mutSummaryMafHypermutators, suffix=doubleSuffix)
    topNDoublesDict = dict([(topHypermutatorGenesArr[:nDoublesToPlot][i][0].split('_')[0], i) for i in range(nDoublesToPlot)]) #complicated line that maps the first substring of the tuple to their index in a sorted list
    genesToPlot = topNDoublesDict.keys()
    
    summaryInfoHypermutatorMaf = process_mut_summary_maf_for_averages(mutSummaryMafHypermutators, genesToPlot, labelMode='hypermutator', suffix=doubleSuffix)
    summaryInfoNonHypermutatorMaf = process_mut_summary_maf_for_averages(mutSummaryMafNonHypermutators, genesToPlot, labelMode='NOT-hypermutator', suffix=doubleSuffix)
    allListOfDicts = summaryInfoHypermutatorMaf + summaryInfoNonHypermutatorMaf
    df = pd.DataFrame(allListOfDicts)
    df['orderingVal'] = df.apply(lambda row: assign_ordering_val(row, topNDoublesDict), axis=1)
    df['isHypermutation'] = df['label'].apply(lambda x: False if 'NOT-' in x else True)
    
    return df

def get_n_multiplets_for_gene(geneDf, caseCol='case'):
    cntrObj = Counter(geneDf[caseCol])
    multipleIds = [i for i in dict(cntrObj).items() if i[1] > 1] #this is probably convoluted and slow; fix!
    nMultipletSamples = len(set(multipleIds))
    return nMultipletSamples

#TODO MOVE THIS ELSEWHERE WHEN BACK AT WORK
def compare_multiplet_ratio_fractions_hyper_vs_non_hyper(mafDf):
    listOfDicts = []
    genes = set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KMT2B', 'KMT5A', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'KMT2A', 'KMT2B', 'KMT2C', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51C', 'RAD51L1', 'RAD51L3', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2'])
    oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
    
    nHypermutatedCases = len(set(mafDf[mafDf['isHypermutator'] == True]['Tumor_Sample_Barcode']))
    nNotHyperCases = len(set(mafDf[mafDf['isHypermutator'] == False]['Tumor_Sample_Barcode']))
        
    for gene in genes:
        geneMaf = mafDf[mafDf['Hugo_Symbol'] == gene]  
        hypermutatorMaf = geneMaf[geneMaf['isHypermutator'] == True]
        nonHypermutatorMaf = geneMaf[geneMaf['isHypermutator'] == False]
        
        hypermuatorOncogenic = hypermutatorMaf[hypermutatorMaf['oncogenic'].isin(oncogenicMutColNames)]
        nonHypermuatorOncogenic = nonHypermutatorMaf[nonHypermutatorMaf['oncogenic'].isin(oncogenicMutColNames)]
        
        listOfDicts.append({
                'gene': gene,
                'nCasesMutated_NonHypermutator': len(set(nonHypermutatorMaf['Tumor_Sample_Barcode'])),
                'nCasesMutated_Hypermutator': len(set(hypermutatorMaf['Tumor_Sample_Barcode'])),
                'nCasesMutated_oncogenic_NonHypermutator': len(set(nonHypermuatorOncogenic['Tumor_Sample_Barcode'])),
                'nCasesMutated_oncogenic_Hypermutator': len(set(hypermuatorOncogenic['Tumor_Sample_Barcode'])),
                'nMultiplet_NonHypermutator': get_n_multiplets_for_gene(nonHypermutatorMaf, caseCol='Tumor_Sample_Barcode'),
                'nMultiplet_Hypermutator': get_n_multiplets_for_gene(hypermutatorMaf, caseCol='Tumor_Sample_Barcode'),
                'nMultiplet_oncogenic_NonHypermutator': get_n_multiplets_for_gene(nonHypermuatorOncogenic, caseCol='Tumor_Sample_Barcode'),
                'nMultiplet_oncogenic_Hypermutator': get_n_multiplets_for_gene(hypermuatorOncogenic, caseCol='Tumor_Sample_Barcode')
                })
    df = pd.DataFrame(listOfDicts)
    
    df['frac_cohort_multiplet_hypermutator'] = df['nMultiplet_Hypermutator'].apply(lambda x: (1.0*x)/nHypermutatedCases)
    df['frac_cohort_multiplet_notHypermutator'] = df['nMultiplet_NonHypermutator'].apply(lambda x: (1.0*x)/nNotHyperCases)
    df['frac_cohort_multiplet_oncogenic_hypermutator'] = df['nMultiplet_oncogenic_Hypermutator'].apply(lambda x: (1.0*x)/nHypermutatedCases)
    df['frac_cohort_multiplet_oncogenic_notHypermutator'] = df['nMultiplet_oncogenic_NonHypermutator'].apply(lambda x: (1.0*x)/nNotHyperCases)
    
    df['oncogenic_multipletRatioHypermutator'] = df.apply(lambda row: 1.0*row['nMultiplet_oncogenic_Hypermutator']/row['nCasesMutated_oncogenic_Hypermutator'] if row['nCasesMutated_oncogenic_Hypermutator'] > 0  else None, axis =1)
    df['oncogenic_multipletRatioNonHypermutator'] = df.apply(lambda row: 1.0*row['nMultiplet_oncogenic_NonHypermutator']/row['nCasesMutated_oncogenic_NonHypermutator'] if row['nCasesMutated_oncogenic_NonHypermutator'] > 0  else None, axis =1)
    df['all_multipletRatioHypermutator'] = df.apply(lambda row: 1.0*row['nMultiplet_Hypermutator']/row['nCasesMutated_Hypermutator'] if row['nCasesMutated_oncogenic_Hypermutator'] > 0  else None, axis =1)
    df['all_multipletRatioNonHypermutator'] = df.apply(lambda row: 1.0*row['nMultiplet_NonHypermutator']/row['nCasesMutated_NonHypermutator'] if row['nCasesMutated_oncogenic_NonHypermutator'] > 0  else None, axis =1)

    return df

def get_glioma_hyper_mmr_gene_fracs(gliomaDf, cases):
    mmrGenes = set(['PMS2', 'MLH1', 'MSH2', 'MSH6'])
    nMMRDoubles = 0
    nMMRMutatedCases = 0
    for case in cases:
        caseMaf = gliomaDf[gliomaDf['Tumor_Sample_Barcode'] == case]
        for gene in mmrGenes:
            geneMaf = caseMaf[caseMaf['Hugo_Symbol'] == gene]
            if geneMaf.shape[0] == 1:
                nMMRMutatedCases += 1
            elif geneMaf.shape[0] > 1:
                nMMRDoubles += 1
                nMMRMutatedCases += 1
                break
    return 1.0*nMMRDoubles/len(hypermutatorCases), 1.0*nMMRDoubles/nMMRMutatedCases
        
def enumerate_recurrently_mutated_tumor_supressors_and_oncogenes(cohortMaf, thresh=.1):
    tumorSupressors = set(['ERRFI1', 'ASXL2', 'PMAIP1', 'ACTG1', 'SUFU', 'FBXO11', 'MEN1', 'FAM58A', 'B2M', 'RB1', 'DUSP22', 'SESN1', 'GPS2', 'RAD51D', 'SMG1', 'CDC73', 'MAP3K1', 'SMARCB1', 'INPP4B', 'PARK2', 'SMAD4', 'CBFB', 'CDH1', 'PPP6C', 'SETDB1', 'SETDB2', 'NF2', 'CDKN2B', 'CDKN2C', 'CDKN2A', 'DDX3X', 'PIK3R1', 'BARD1', 'PDS5B', 'KLF4', 'SPRED1', 'VHL', 'SMAD2', 'PMS1', 'PMS2', 'SETD2', 'GATA3', 'TBL1XR1', 'MUTYH', 'SOCS1', 'FAM175A', 'ROBO1', 'ARID1B', 'ARID1A', 'TCF7L2', 'STK11', 'FOXA1', 'PTEN', 'FAT1', 'FAS', 'CYLD', 'MAX', 'SH2D1A', 'APC', 'NTHL1', 'CTCF', 'KDM5C', 'KMT2C', 'ZFHX3', 'FOXP1', 'PIGA', 'CDKN1B', 'CDKN1A', 'FUBP1', 'MSH2', 'ID3', 'TNFRSF14', 'TRAF3', 'EP400', 'BRIP1', 'ARID4A', 'ARID4B', 'XRCC2', 'DAXX', 'SDHAF2', 'ASXL1', 'AMER1', 'RASA1', 'EGR1', 'MST1', 'SOX17', 'RUNX1', 'PIK3R3', 'NCOR1', 'NF1', 'JAK1', 'PTPRD', 'CHEK2', 'CHEK1', 'SMC1A', 'TMEM127', 'STAG1', 'RAD51', 'TCF3', 'STAG2', 'ARID2', 'RAD50', 'RNF43', 'PARP1', 'BLM', 'CUX1', 'RECQL', 'RAD21', 'PTPN2', 'PTPN1', 'SLX4', 'INHA', 'PAX5', 'IRF1', 'TP53', 'HLA-A', 'IRF8', 'CBL', 'TOP1', 'SHQ1', 'PRDM1', 'NSD1', 'ATXN2', 'CREBBP', 'HDAC4', 'SESN2', 'PPP2R1A', 'EPHA7', 'ATM', 'EPHA3', 'POT1', 'SMAD3', 'MOB3B', 'TBX3', 'POLE', 'ATR', 'FANCD2', 'FH', 'BCORL1', 'SOX9', 'IKZF3', 'TSC1', 'TP63', 'MRE11A', 'SDHC', 'BTG1', 'POLD1', 'CIITA', 'SMC3', 'SAMHD1', 'RTEL1', 'ECT2L', 'PIK3R2', 'CRBN', 'FANCC', 'NBN', 'FANCA', 'HLA-B', 'RECQL4', 'DUSP4', 'ERCC2', 'FBXW7', 'TGFBR2', 'TGFBR1', 'MSH3', 'RBM15', 'TET1', 'TET3', 'SESN3', 'MGA', 'LTB', 'FOXL2', 'SH2B3', 'BCOR', 'HIST1H1D', 'ATRX', 'EP300', 'RAD51C', 'RAD51B', 'HIST1H1B', 'TNFAIP3', 'DICER1', 'ARID5B', 'LATS2', 'FOXO1', 'KEAP1', 'EZH2', 'SP140', 'NKX3-1', 'PBRM1', 'PALB2', 'CIC', 'BRCA1', 'DTX1', 'FLCN', 'SPEN', 'CD58', 'ERCC3', 'ERCC4', 'MSH6', 'BCL11B', 'BMPR1A', 'ERF', 'BRCA2', 'NOTCH2', 'EED', 'MITF', 'ELF3', 'SMARCA4', 'BBC3', 'ANKRD11', 'CEBPA', 'BCL2L11', 'AXIN2', 'AXIN1', 'CDK12', 'ESCO2', 'MLH1', 'SDHB', 'MED12', 'HNF1A', 'RYBP', 'ATP6V1B2', 'DNMT3B', 'KMT2B', 'KMT2A', 'DNMT3A', 'NFKBIA', 'TRAF5', 'KMT2D', 'SPOP', 'RBM10', 'P2RY8', 'TP53BP1', 'TSC2', 'KDM6A', 'EPCAM', 'PHOX2B', 'NPM1', 'BCL10', 'LATS1', 'HOXB13', 'ARID3A', 'PTPRT', 'PTPRS', 'INPPL1', 'NOTCH4', 'TET2', 'NOTCH1', 'CASP8', 'NOTCH3', 'GRIN2A', 'MAP2K4', 'WT1', 'BACH2', 'SDHA', 'BAP1', 'PTCH1', 'SDHD'])
    occurenceCounter, rankingDict, fractionalDict = maf_analysis_utils.enumerate_top_n_oncogenic_mutated_genes_across_cohort(cohortMaf, n=50)
    recurrentTumorSupressors = []
    recurrentOncogenes = []
    for key, value in fractionalDict.items():
        if value > .1:
            if key in tumorSupressors:
                recurrentTumorSupressors.append(key)
            else:
                recurrentOncogenes.append(key) 
    return set(recurrentTumorSupressors), set(recurrentOncogenes)

def pan_hypermutated_cancer_type_double_prevalence_analysis(maf, cancerTypes = ['Endometrial Cancer', 'Glioma', 'Colorectal Cancer', 'Bladder Cancer', 'Melanoma', 'Non-Small Cell Lung Cancer']):
    tumorSupressors = set(['ERRFI1', 'ASXL2', 'PMAIP1', 'ACTG1', 'SUFU', 'FBXO11', 'MEN1', 'FAM58A', 'B2M', 'RB1', 'DUSP22', 'SESN1', 'GPS2', 'RAD51D', 'SMG1', 'CDC73', 'MAP3K1', 'SMARCB1', 'INPP4B', 'PARK2', 'SMAD4', 'CBFB', 'CDH1', 'PPP6C', 'SETDB1', 'SETDB2', 'NF2', 'CDKN2B', 'CDKN2C', 'CDKN2A', 'DDX3X', 'PIK3R1', 'BARD1', 'PDS5B', 'KLF4', 'SPRED1', 'VHL', 'SMAD2', 'PMS1', 'PMS2', 'SETD2', 'GATA3', 'TBL1XR1', 'MUTYH', 'SOCS1', 'FAM175A', 'ROBO1', 'ARID1B', 'ARID1A', 'TCF7L2', 'STK11', 'FOXA1', 'PTEN', 'FAT1', 'FAS', 'CYLD', 'MAX', 'SH2D1A', 'APC', 'NTHL1', 'CTCF', 'KDM5C', 'KMT2C', 'ZFHX3', 'FOXP1', 'PIGA', 'CDKN1B', 'CDKN1A', 'FUBP1', 'MSH2', 'ID3', 'TNFRSF14', 'TRAF3', 'EP400', 'BRIP1', 'ARID4A', 'ARID4B', 'XRCC2', 'DAXX', 'SDHAF2', 'ASXL1', 'AMER1', 'RASA1', 'EGR1', 'MST1', 'SOX17', 'RUNX1', 'PIK3R3', 'NCOR1', 'NF1', 'JAK1', 'PTPRD', 'CHEK2', 'CHEK1', 'SMC1A', 'TMEM127', 'STAG1', 'RAD51', 'TCF3', 'STAG2', 'ARID2', 'RAD50', 'RNF43', 'PARP1', 'BLM', 'CUX1', 'RECQL', 'RAD21', 'PTPN2', 'PTPN1', 'SLX4', 'INHA', 'PAX5', 'IRF1', 'TP53', 'HLA-A', 'IRF8', 'CBL', 'TOP1', 'SHQ1', 'PRDM1', 'NSD1', 'ATXN2', 'CREBBP', 'HDAC4', 'SESN2', 'PPP2R1A', 'EPHA7', 'ATM', 'EPHA3', 'POT1', 'SMAD3', 'MOB3B', 'TBX3', 'POLE', 'ATR', 'FANCD2', 'FH', 'BCORL1', 'SOX9', 'IKZF3', 'TSC1', 'TP63', 'MRE11A', 'SDHC', 'BTG1', 'POLD1', 'CIITA', 'SMC3', 'SAMHD1', 'RTEL1', 'ECT2L', 'PIK3R2', 'CRBN', 'FANCC', 'NBN', 'FANCA', 'HLA-B', 'RECQL4', 'DUSP4', 'ERCC2', 'FBXW7', 'TGFBR2', 'TGFBR1', 'MSH3', 'RBM15', 'TET1', 'TET3', 'SESN3', 'MGA', 'LTB', 'FOXL2', 'SH2B3', 'BCOR', 'HIST1H1D', 'ATRX', 'EP300', 'RAD51C', 'RAD51B', 'HIST1H1B', 'TNFAIP3', 'DICER1', 'ARID5B', 'LATS2', 'FOXO1', 'KEAP1', 'EZH2', 'SP140', 'NKX3-1', 'PBRM1', 'PALB2', 'CIC', 'BRCA1', 'DTX1', 'FLCN', 'SPEN', 'CD58', 'ERCC3', 'ERCC4', 'MSH6', 'BCL11B', 'BMPR1A', 'ERF', 'BRCA2', 'NOTCH2', 'EED', 'MITF', 'ELF3', 'SMARCA4', 'BBC3', 'ANKRD11', 'CEBPA', 'BCL2L11', 'AXIN2', 'AXIN1', 'CDK12', 'ESCO2', 'MLH1', 'SDHB', 'MED12', 'HNF1A', 'RYBP', 'ATP6V1B2', 'DNMT3B', 'KMT2B', 'KMT2A', 'DNMT3A', 'NFKBIA', 'TRAF5', 'KMT2D', 'SPOP', 'RBM10', 'P2RY8', 'TP53BP1', 'TSC2', 'KDM6A', 'EPCAM', 'PHOX2B', 'NPM1', 'BCL10', 'LATS1', 'HOXB13', 'ARID3A', 'PTPRT', 'PTPRS', 'INPPL1', 'NOTCH4', 'TET2', 'NOTCH1', 'CASP8', 'NOTCH3', 'GRIN2A', 'MAP2K4', 'WT1', 'BACH2', 'SDHA', 'BAP1', 'PTCH1', 'SDHD'])
    mutSummaryMafs = []
    for cType in cancerTypes:
        print cType
        cTypeMaf = maf[maf['cancer_type'] == cType]
        recurrentTumorSupressors, recurrentOncogenes = enumerate_recurrently_mutated_tumor_supressors_and_oncogenes(cTypeMaf)
        mutSummaryMaf = compare_multiplet_ratio_fractions_hyper_vs_non_hyper(cTypeMaf)
        mutSummaryMaf['geneClassification'] = mutSummaryMaf['gene'].apply(lambda x: cType + '_RecurrentTumorSupressor' if x in recurrentTumorSupressors
                     else cType + '_RecurrentOncogene' if x in recurrentOncogenes
                     else cType + '_NonRecurrentTumorSupressor' if x in tumorSupressors
                     else cType + '_NonRecurrentOncogene')
        mutSummaryMaf['cancerType'] = cType
        mutSummaryMafs.append(mutSummaryMaf)
    return mutSummaryMafs

def enumerate_double_mutated_genes_and_allele_prevalence(mafDf, oncogenicOnly=True):
    
    mafDf = mafDf[mafDf['HGVSp_Short'].notnull()]
    listOfDicts = []
    genes = set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KMT2B', 'KMT5A', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'KMT2A', 'KMT2B', 'KMT2C', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51C', 'RAD51L1', 'RAD51L3', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2'])
    
    oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
    if oncogenicOnly: #if we only do the 
        mafDf = mafDf[mafDf['oncogenic'].isin(oncogenicMutColNames)]

    genes = set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KMT2B', 'KMT5A', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'KMT2A', 'KMT2B', 'KMT2C', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51C', 'RAD51L1', 'RAD51L3', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2'])
    cntr = 0
    
    
    for gene in genes:
        if cntr % 25==0: print cntr
        cntr += 1
        geneMaf = mafDf[mafDf['Hugo_Symbol'] == gene]  
        
        for case in set(geneMaf['Tumor_Sample_Barcode']):
            caseMaf = geneMaf[geneMaf['Tumor_Sample_Barcode'] == case]
            if caseMaf.shape[0] > 1:
                for index, row in caseMaf.iterrows():
                    allele = row['HGVSp_Short']
                    indel = False
                    if '*' in allele: indel = True
                    listOfDicts.append(
                            {'Tumor_Sample_Barcode': case,
                             'Gene': gene,
                             'allele': allele,
                             'indel': indel}
                            )
    
    return pd.DataFrame(listOfDicts)
      
    
    
def filter_and_prep_maf_for_analyis(mafDf, cType):
    
    cHyperIds = analysis_utils.get_ids_by_hypermutant_status(hypermutantIdDir=pathPrefix + '/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/hypermutationStatusIds', cancerType=cType, hypermutantStatus = 'Hypermutated')
    cNormalIds = analysis_utils.get_ids_by_hypermutant_status(hypermutantIdDir=pathPrefix +'/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/hypermutationStatusIds', cancerType=cType, hypermutantStatus = 'Normal')
    cMaf = mafDf[mafDf['Tumor_Sample_Barcode'].isin(cHyperIds | cNormalIds)]
    cMaf['sampleMutStatus'] = cMaf['Tumor_Sample_Barcode'].apply(lambda x:
        'hypermutated' if x in cHyperIds
        else 'normal'
        )
    
    return cMaf
    
    

maf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/hypermutationAnalysisProjectData/hypermutationAnalysisProjectMaf.tsv')
sigs = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/hypermutationAnalysisProjectData/hypermutationAnalysisProjectSignatures.tsv')

#FIX MLL IN THE MAF

#TODO FOR NEXT WEEK
#MARK THE DOUBLETS AHEAD OF TIME THEN GO TO Them and do processing and analysis etc


maf = analysis_utils.fix_mll_genes(maf)
dfInfo = filter_and_prep_maf_for_analyis(maf, 'Endometrial Cancer')
dfSummary = enumerate_double_mutated_genes_and_allele_prevalence(dfInfo, oncogenicOnly=True)
alleleSummaryDict = dict(Counter(dfSummary['allele']))
geneSummaryDict = dict(Counter(dfSummary['Gene']))
alleleRecurrenceDisplayThresh = 10
dfSummary['AlleleCount'] = dfSummary['allele'].apply(lambda x: alleleSummaryDict[x])
dfSummary['label'] = dfSummary['allele'].apply(lambda x: x if alleleSummaryDict[x] > alleleRecurrenceDisplayThresh else None)
dfSummary['orderingVal'] = dfSummary['Gene'].apply(lambda x: geneSummaryDict[x])
dfSummary = dfSummary[dfSummary['orderingVal'] > 15]
dfSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/endometrialDoubleMutationPlotting.tsv', index=False, sep='\t')

print Counter(dfSummary['allele']).most_common(10)

maf.columns.values











endometrialMaf = mafSansIntrons[mafSansIntrons['cancer_type'] == 'Endometrial Cancer']
mutSummaryMafEndometrial = compare_multiplet_ratio_fractions_hyper_vs_non_hyper(endometrialMaf)

gliomaMaf = mafSansIntrons[mafSansIntrons['cancer_type'] == 'Glioma']
mutSummaryMafGlioma = compare_multiplet_ratio_fractions_hyper_vs_non_hyper(gliomaMaf)

colorectalMaf = mafSansIntrons[mafSansIntrons['cancer_type'] == 'Colorectal Cancer']
mutSummaryMafColorectal = compare_multiplet_ratio_fractions_hyper_vs_non_hyper(colorectalMaf)

#WE USE THE FILTERED MAFs for this analysis because we dont want to deal with spurious double oncogenic mutations


#OLDER CODE
#DO SPECIFIC STUFF TO SET WHICH LABELS WE DO AND DONT DISPLAY
mutSummaryMafEndometrial['textLabel'] = mutSummaryMafEndometrial.apply(lambda row: row['gene'] if ((row['oncogenic_multipletRatioHypermutator'] > .25) & (row['frac_cohort_multiplet_oncogenic_hypermutator'] > .1)) else None, axis=1)
mutSummaryMafEndometrial.to_csv('~/Desktop/WORK/dataForLocalPlotting/endometrial_double_summary.tsv', index=False, sep='\t')

mutSummaryMafColorectal['textLabel'] = mutSummaryMafColorectal.apply(lambda row: row['gene'] if ((row['oncogenic_multipletRatioHypermutator'] > .25) & (row['frac_cohort_multiplet_oncogenic_hypermutator'] > .1)) else None, axis=1)
mutSummaryMafColorectal.to_csv('~/Desktop/WORK/dataForLocalPlotting/colorectal_double_summary.tsv', index=False, sep='\t')

hypermutatorCases = set(gliomaMaf[gliomaMaf['isHypermutator'] == True]['Tumor_Sample_Barcode'])
notHypermutatorCases = set(gliomaMaf[gliomaMaf['isHypermutator'] == False]['Tumor_Sample_Barcode'])
f1, f2 = get_glioma_hyper_mmr_gene_fracs(gliomaMaf, hypermutatorCases)
n1, n2 = get_glioma_hyper_mmr_gene_fracs(gliomaMaf, notHypermutatorCases)
gliomaGeneFracDf = pd.DataFrame([{'frac_cohort_multiplet_hypermutator': f1, 'frac_cohort_multiplet_notHypermutator':n1, 'all_multipletRatioHypermutator':f2, 'all_multipletRatioNonHypermutator':n2, 'gene': 'ALL_MMR_GENES'}])
mutSummaryMafGlioma = pd.concat([mutSummaryMafGlioma, gliomaGeneFracDf])

mutSummaryMafGlioma['textLabel'] = mutSummaryMafGlioma.apply(lambda row: row['gene'] if ((row['all_multipletRatioHypermutator'] > .5) & (row['frac_cohort_multiplet_hypermutator'] > .1)) else None, axis=1)
mutSummaryMafGlioma.to_csv('~/Desktop/WORK/dataForLocalPlotting/glioma_double_summary.tsv', index=False, sep='\t')

mutSummaryMafGlioma[mutSummaryMafGlioma['gene'] == 'MSH6']['frac_cohort_multiplet_notHypermutator']

tumorSupressors = set(['ERRFI1', 'ASXL2', 'PMAIP1', 'ACTG1', 'SUFU', 'FBXO11', 'MEN1', 'FAM58A', 'B2M', 'RB1', 'DUSP22', 'SESN1', 'GPS2', 'RAD51D', 'SMG1', 'CDC73', 'MAP3K1', 'SMARCB1', 'INPP4B', 'PARK2', 'SMAD4', 'CBFB', 'CDH1', 'PPP6C', 'SETDB1', 'SETDB2', 'NF2', 'CDKN2B', 'CDKN2C', 'CDKN2A', 'DDX3X', 'PIK3R1', 'BARD1', 'PDS5B', 'KLF4', 'SPRED1', 'VHL', 'SMAD2', 'PMS1', 'PMS2', 'SETD2', 'GATA3', 'TBL1XR1', 'MUTYH', 'SOCS1', 'FAM175A', 'ROBO1', 'ARID1B', 'ARID1A', 'TCF7L2', 'STK11', 'FOXA1', 'PTEN', 'FAT1', 'FAS', 'CYLD', 'MAX', 'SH2D1A', 'APC', 'NTHL1', 'CTCF', 'KDM5C', 'KMT2C', 'ZFHX3', 'FOXP1', 'PIGA', 'CDKN1B', 'CDKN1A', 'FUBP1', 'MSH2', 'ID3', 'TNFRSF14', 'TRAF3', 'EP400', 'BRIP1', 'ARID4A', 'ARID4B', 'XRCC2', 'DAXX', 'SDHAF2', 'ASXL1', 'AMER1', 'RASA1', 'EGR1', 'MST1', 'SOX17', 'RUNX1', 'PIK3R3', 'NCOR1', 'NF1', 'JAK1', 'PTPRD', 'CHEK2', 'CHEK1', 'SMC1A', 'TMEM127', 'STAG1', 'RAD51', 'TCF3', 'STAG2', 'ARID2', 'RAD50', 'RNF43', 'PARP1', 'BLM', 'CUX1', 'RECQL', 'RAD21', 'PTPN2', 'PTPN1', 'SLX4', 'INHA', 'PAX5', 'IRF1', 'TP53', 'HLA-A', 'IRF8', 'CBL', 'TOP1', 'SHQ1', 'PRDM1', 'NSD1', 'ATXN2', 'CREBBP', 'HDAC4', 'SESN2', 'PPP2R1A', 'EPHA7', 'ATM', 'EPHA3', 'POT1', 'SMAD3', 'MOB3B', 'TBX3', 'POLE', 'ATR', 'FANCD2', 'FH', 'BCORL1', 'SOX9', 'IKZF3', 'TSC1', 'TP63', 'MRE11A', 'SDHC', 'BTG1', 'POLD1', 'CIITA', 'SMC3', 'SAMHD1', 'RTEL1', 'ECT2L', 'PIK3R2', 'CRBN', 'FANCC', 'NBN', 'FANCA', 'HLA-B', 'RECQL4', 'DUSP4', 'ERCC2', 'FBXW7', 'TGFBR2', 'TGFBR1', 'MSH3', 'RBM15', 'TET1', 'TET3', 'SESN3', 'MGA', 'LTB', 'FOXL2', 'SH2B3', 'BCOR', 'HIST1H1D', 'ATRX', 'EP300', 'RAD51C', 'RAD51B', 'HIST1H1B', 'TNFAIP3', 'DICER1', 'ARID5B', 'LATS2', 'FOXO1', 'KEAP1', 'EZH2', 'SP140', 'NKX3-1', 'PBRM1', 'PALB2', 'CIC', 'BRCA1', 'DTX1', 'FLCN', 'SPEN', 'CD58', 'ERCC3', 'ERCC4', 'MSH6', 'BCL11B', 'BMPR1A', 'ERF', 'BRCA2', 'NOTCH2', 'EED', 'MITF', 'ELF3', 'SMARCA4', 'BBC3', 'ANKRD11', 'CEBPA', 'BCL2L11', 'AXIN2', 'AXIN1', 'CDK12', 'ESCO2', 'MLH1', 'SDHB', 'MED12', 'HNF1A', 'RYBP', 'ATP6V1B2', 'DNMT3B', 'KMT2B', 'KMT2A', 'DNMT3A', 'NFKBIA', 'TRAF5', 'KMT2D', 'SPOP', 'RBM10', 'P2RY8', 'TP53BP1', 'TSC2', 'KDM6A', 'EPCAM', 'PHOX2B', 'NPM1', 'BCL10', 'LATS1', 'HOXB13', 'ARID3A', 'PTPRT', 'PTPRS', 'INPPL1', 'NOTCH4', 'TET2', 'NOTCH1', 'CASP8', 'NOTCH3', 'GRIN2A', 'MAP2K4', 'WT1', 'BACH2', 'SDHA', 'BAP1', 'PTCH1', 'SDHD'])



v = pan_hypermutated_cancer_type_double_prevalence_analysis(mafSansIntrons)

df = pd.concat(v)
df['displayGeneName'] = df.apply(lambda row: row['gene'] if row['frac_cohort_multiplet_oncogenic_hypermutator'] > .1 
  else row['gene'] if 'Non' not in row['geneClassification']
  else None, axis=1)

orderingDict1 = {'Endometrial Cancer': 0, 'Colorectal Cancer': 1, 'Glioma': 2, 'Melanoma': 3, 'Bladder Cancer': 4, 'Non-Small Cell Lung Cancer':5}
orderingDict2 = {'RecurrentTumorSupressor': .1, 'RecurrentOncogene': .2, 'NonRecurrentTumorSupressor': .3, 'NonRecurrentOncogene':.4}
df['orderingVal'] = df['geneClassification'].apply(lambda x:
    orderingDict1[x.split('_')[0]] + orderingDict2[x.split('_')[1]] 
    )

df['geneClassificationPanCan'] = df['geneClassification'].apply(lambda x: 
  'NonRecurrentTumorSupressor' if 'NonRecurrentTumorSupressor' in x
  else 'NonRecurrentOncogene' if 'NonRecurrentOncogene' in x
  else 'RecurrentTumorSupressor' if 'RecurrentTumorSupressor' in x
  else 'RecurrentOncogene' if 'RecurrentOncogene' in x 
  else None
   )

df['geneAndCohort'] = df.apply(lambda row: row['gene'] + '_' + row['cancerType'], axis=1)
df['geneAndCohortDisplay'] = df.apply(lambda row: row['geneAndCohort'] if row['frac_cohort_multiplet_oncogenic_hypermutator'] > .18 else None, axis=1)
df['displayGlioma'] = df.apply(lambda row: row['gene'] if row['gene'] == 'MSH6' else row['gene'] if row['gene'] == 'MSH2' or row['gene'] == 'PMS2' or row['gene'] == 'MLH1' else row['gene'] if row['gene'] == 'TP53' or row['gene'] == 'NF1' or row['gene'] == 'ATRX' else  None, axis=1)

df.to_csv('~/Desktop/WORK/dataForLocalPlotting/pancanDoubleSummary.tsv', index=False, sep='\t')


df['geneClassification']


print df[(df['cancerType'] == 'Glioma') & (df['gene'] == 'TERT')]['nMultiplet_Hypermutator']














#########ENDOMETRIAL CANCER
endometrialMaf = mafSansIntrons[mafSansIntrons['cancer_type'] == 'Endometrial Cancer']
mutSummaryMafEndometrial = double_mutation_analysis_util.create_double_mutation_summary_vanilla_maf(endometrialMaf)

endometrialDoubleSummary = create_hypermutation_prevalence_summary(mutSummaryMafEndometrial, nDoublesToPlot = 10, doubleSuffix='_double')
endometrialDoubleWithOncogenicSummary =create_hypermutation_prevalence_summary(mutSummaryMafEndometrial, nDoublesToPlot = 10, doubleSuffix='_doubleWithOncogenic')
endometrialDoubleOncogenicSummary = create_hypermutation_prevalence_summary(mutSummaryMafEndometrial, nDoublesToPlot = 10, doubleSuffix='_doubleOncogenic')

endometrialDoubleSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/endometrial_double_mutation_incidence.tsv', sep='\t', index=False)
endometrialDoubleWithOncogenicSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/endometrial_double_with_oncogenic_mutation_incidence.tsv', sep='\t', index=False)
endometrialDoubleOncogenicSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/endometrial_double_with_both_oncogenic_mutation_incidence.tsv', sep='\t', index=False)

##########COLORECTAL CANCER
colorectalMaf = mafSansIntrons[mafSansIntrons['cancer_type'] == 'Colorectal Cancer']
mutSummaryMafColorectal = double_mutation_analysis_util.create_double_mutation_summary_vanilla_maf(colorectalMaf)

colorectalDoubleSummary = create_hypermutation_prevalence_summary(mutSummaryMafColorectal, nDoublesToPlot = 10, doubleSuffix='_double')
colorectalDoubleWithOncogenicSummary =create_hypermutation_prevalence_summary(mutSummaryMafColorectal, nDoublesToPlot = 10, doubleSuffix='_doubleWithOncogenic')
colorectalDoubleOncogenicSummary = create_hypermutation_prevalence_summary(mutSummaryMafColorectal, nDoublesToPlot = 10, doubleSuffix='_doubleOncogenic')

colorectalDoubleSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/colorectal_double_mutation_incidence.tsv', sep='\t', index=False)
colorectalDoubleWithOncogenicSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/colorectal_double_with_oncogenic_mutation_incidence.tsv', sep='\t', index=False)
colorectalDoubleOncogenicSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/colorectal_double_with_both_oncogenic_mutation_incidence.tsv', sep='\t', index=False)

##########GLIOMA
gliomaMaf = mafSansIntrons[mafSansIntrons['cancer_type'] == 'Glioma']
mutSummaryMafGlioma = double_mutation_analysis_util.create_double_mutation_summary_vanilla_maf(gliomaMaf)

gliomaDoubleSummary = create_hypermutation_prevalence_summary(mutSummaryMafGlioma, nDoublesToPlot = 10, doubleSuffix='_double')
gliomaDoubleWithOncogenicSummary =create_hypermutation_prevalence_summary(mutSummaryMafGlioma, nDoublesToPlot = 10, doubleSuffix='_doubleWithOncogenic')
gliomaDoubleOncogenicSummary = create_hypermutation_prevalence_summary(mutSummaryMafGlioma, nDoublesToPlot = 10, doubleSuffix='_doubleOncogenic')

gliomaDoubleSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/glioma_double_mutation_incidence.tsv', sep='\t', index=False)
gliomaDoubleWithOncogenicSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/glioma_double_with_oncogenic_mutation_incidence.tsv', sep='\t', index=False)
gliomaDoubleOncogenicSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/glioma_double_with_both_oncogenic_mutation_incidence.tsv', sep='\t', index=False)

###############BLADDER CANCER

bladderMaf = mafSansIntrons[mafSansIntrons['cancer_type'] == 'Bladder Cancer']
mutSummaryMafBladder = double_mutation_analysis_util.create_double_mutation_summary_vanilla_maf(bladderMaf)

bladderDoubleSummary = create_hypermutation_prevalence_summary(mutSummaryMafBladder, nDoublesToPlot = 10, doubleSuffix='_double')
bladderDoubleWithOncogenicSummary =create_hypermutation_prevalence_summary(mutSummaryMafBladder, nDoublesToPlot = 10, doubleSuffix='_doubleWithOncogenic')
bladderDoubleOncogenicSummary = create_hypermutation_prevalence_summary(mutSummaryMafBladder, nDoublesToPlot = 10, doubleSuffix='_doubleOncogenic')

bladderDoubleSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/bladder_double_mutation_incidence.tsv', sep='\t', index=False)
bladderDoubleWithOncogenicSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/bladder_double_with_oncogenic_mutation_incidence.tsv', sep='\t', index=False)
bladderDoubleOncogenicSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/bladder_double_with_both_oncogenic_mutation_incidence.tsv', sep='\t', index=False)

################MELANOMA

melanomaMaf = mafSansIntrons[mafSansIntrons['cancer_type'] == 'Melanoma']
mutSummaryMafMelanoma = double_mutation_analysis_util.create_double_mutation_summary_vanilla_maf(melanomaMaf)

melanomaDoubleSummary = create_hypermutation_prevalence_summary(mutSummaryMafMelanoma, nDoublesToPlot = 10, doubleSuffix='_double')
melanomaDoubleWithOncogenicSummary =create_hypermutation_prevalence_summary(mutSummaryMafMelanoma, nDoublesToPlot = 10, doubleSuffix='_doubleWithOncogenic')
melanomaDoubleOncogenicSummary = create_hypermutation_prevalence_summary(mutSummaryMafMelanoma, nDoublesToPlot = 10, doubleSuffix='_doubleOncogenic')

melanomaDoubleSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/melanoma_double_mutation_incidence.tsv', sep='\t', index=False)
melanomaDoubleWithOncogenicSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/melanoma_double_with_oncogenic_mutation_incidence.tsv', sep='\t', index=False)
melanomaDoubleOncogenicSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/melanoma_double_with_both_oncogenic_mutation_incidence.tsv', sep='\t', index=False)






