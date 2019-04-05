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
import math
import mutation_modeling_util

def convert_dict_to_df(d):
    genes = d.keys()
    listOfDicts = []
    for gene in genes:
        localD = d[gene]
        for quadNuc in localD.keys():
            listOfDicts.append({'Hugo_Symbol': gene, 'quadNuc': quadNuc,
                                'totalNmut': localD[quadNuc]['totalNmut'], 'nOncogenicMut': localD[quadNuc]['nOncogenicMut']})
    return pd.DataFrame(listOfDicts)

#CALCULATES THE DIFFERENTIAL EXPECTED oncogenic mut burden as if each signature were dominant alone
def calculate_oncogenic_mut_susceptibility_of_genes_by_signature(oncogenicSDict):
    listOfDicts = []
    sigNames = ['Signature.' + str(i) for i in range(1,31)]
    for i in range(1,31):
        curSig = 'Signature.' + str(i)
        d = {}
        for s in sigNames:
            d[s] = 0
        d[curSig] = 1
        #PRETEND we got a case with 100% signature i on the decomposition
        quadNucFractions = mutation_modeling_util.get_quadnuc_fracs_given_decomposition(d, spectraPath = pathPrefix + '/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt')
        v = mutation_modeling_util.get_expected_oncogenic_val_given_quadnuc_fractions(quadNucFractions, oncogenicSDict)
        listOfDicts.append({'Signature_Name': curSig, 'ExpectedFracOfMutsOncogenic': v})
    return pd.DataFrame(listOfDicts)



#GIVen a gene, oncogenic info and a signature returns the fraction of mutations that will be an oncogenic mutation in that gene
def calculate_gene_specific_oncogenic_susceptibility_by_signature(signature, oncogenicSDict, gene):
    sigNames = ['Signature.' + str(i) for i in range(1,31)]
    decompositionDict = {}
    for s in sigNames:
        decompositionDict[s] = 0 
    decompositionDict[signature] = 1
    quadNucFractions = mutation_modeling_util.get_quadnuc_fracs_given_decomposition(decompositionDict, spectraPath = pathPrefix + '/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt')
    oncogenicExpected = mutation_modeling_util.get_expected_oncogenic_val_given_quadnuc_fractions_and_gene(quadNucFractions, oncogenicSDict, 'IMPACT_468', gene)
    return oncogenicExpected

def get_top_n_most_mutable_genes_by_signature(signature, oncogenicSDict, n):
    l= []
    impact468genes = set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'AMER1', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'GTF2I', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MLL4', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SETD8', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2']) 
    cntr = 0
    for gene in impact468genes:
        cntr += 1
        if cntr%10 == 0: print cntr
        v = calculate_gene_specific_oncogenic_susceptibility_by_signature(signature, oncogenicSDict, gene)
        l.append((gene, v))
    print 'done calculating top n genes mut by signature'
    return sorted(l, key= lambda tup: tup[1], reverse=True)[:n]
    
def get_genes_mutated_in_cohort_above_threshold(maf, thresh=.1):
    oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
    maf = maf[maf['oncogenic'].isin(oncogenicMutColNames)]
    genesWithEnoughMutFreq = []
    nCases = len(set(maf['Tumor_Sample_Barcode']))
    cntr = 0
    for gene in set(maf['Hugo_Symbol']):
        cntr += 1
        if cntr%10 == 0: print cntr
        nObservedMutInCohort = maf[maf['Hugo_Symbol'] == gene].shape[0]
        if 1.0*nObservedMutInCohort/nCases > thresh:
            genesWithEnoughMutFreq.append(gene)
    print 'done with cohort summarizing'
    return set(genesWithEnoughMutFreq)
        



#TODO OBSERVED VS EXPECTED MUTS IN GENES
def compare_observed_vs_expected_gene_oncogenic_mut_prevalence_across_cohort(cohortMaf, oncogenicSDict, sigs, genes=None):
    panelGeneDict = {
            'IMPACT_368': set(['ABL1', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BAP1', 'BARD1', 'BBC3', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CARD11', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79B', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSF1R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EP300', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERG', 'ESR1', 'ETV1', 'ETV6', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXP1', 'FUBP1', 'GATA1', 'GATA2', 'GATA3', 'GNA11', 'GNAQ', 'GNAS', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3B', 'HNF1A', 'HRAS', 'ICOSLG', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INPP4A', 'INPP4B', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAPK1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH6', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOR1', 'NF1', 'NF2', 'NFE2L2', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTRK1', 'NTRK2', 'NTRK3', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLE', 'PPP2R1A', 'PRDM1', 'PRKAR1A', 'PTCH1', 'PTEN', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAC1', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RUNX1', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SETD2', 'SF3B1', 'SH2D1A', 'SHQ1', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SOCS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SRC', 'STAG2', 'STK11', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TBX3', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP63', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'VHL', 'VTCN1', 'WT1', 'XIAP', 'XPO1', 'YAP1', 'YES1']),
            'IMPACT_410': set(['ABL1', 'ACVR1', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'EP300', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH6', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPM1D', 'PPP2R1A', 'PPP6C', 'PRDM1', 'PRKAR1A', 'PTCH1', 'PTEN', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAD21', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RUNX1', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHQ1', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SOCS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'VEGFA', 'VHL', 'VTCN1', 'WT1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2']),
            'IMPACT_468': set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'AMER1', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'GTF2I', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MLL4', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SETD8', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2']) 
            }
    impactVersionDict = {'IM3': 'IMPACT_368', 'IM5': 'IMPACT_410', 'IM6': 'IMPACT_468'}
    oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
    
    if genes == None:
        genes = panelGeneDict['IMPACT_468']
    
    cohortMafSnpsOnly = cohortMaf[(cohortMaf['Variant_Type'] == 'SNP')]
    cohortMafOncogenicSnpsOnly = cohortMafSnpsOnly[cohortMafSnpsOnly['oncogenic'].isin(oncogenicMutColNames)]
    listOfDicts = []
    cntr = 0
    for gene in genes:
        print 'analyzing gene ', cntr, 'th gene: ', gene
        cntr += 1
        curNumObserved = 0
        curNumExpected = 0
        for case in set(cohortMaf['Tumor_Sample_Barcode']):
            nCaseSnps = cohortMafSnpsOnly[cohortMafSnpsOnly['Tumor_Sample_Barcode'] == case].shape[0]
            casePanel = case[-3:]
            if casePanel not in impactVersionDict:
                pass
            else:
                impactVersion = impactVersionDict[casePanel]
                if gene not in panelGeneDict[impactVersion]: #IF the current case's gene panel didnt include the current gene skip it
                    curNumObserved += 0
                    curNumExpected += 0
                else:
                    nOncogenicSnpsInGeneInCase = cohortMafOncogenicSnpsOnly[(cohortMafOncogenicSnpsOnly['Tumor_Sample_Barcode'] == case) & (cohortMafOncogenicSnpsOnly['Hugo_Symbol'] == gene)].shape[0]
                    curNumObserved += nOncogenicSnpsInGeneInCase
                    
                    #CALCULATE THE N expected oncogenic mutations to a gene expected in a single case
                    caseSigs = sigs[sigs['Tumor_Sample_Barcode'] == case]
                    signatureColumnNames = ['mean_' + str(i) for i in range(1,31)]
                    signatureColsOnly = caseSigs[signatureColumnNames]
                    renameDict = dict([('mean_' + str(i), 'Signature.' + str(i)) for i in range(1,31)])
                    signatureColsOnly = signatureColsOnly.rename(columns=renameDict)
                    
                    decompositionDict = signatureColsOnly.to_dict(orient='records')[0]
                    quadNucFractions = mutation_modeling_util.get_quadnuc_fracs_given_decomposition(decompositionDict, spectraPath = pathPrefix + '/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt')
                    oncogenicExpected = mutation_modeling_util.get_expected_oncogenic_val_given_quadnuc_fractions_and_gene(quadNucFractions, oncogenicSDict, impactVersion, gene)
                    curNumExpected += nCaseSnps*oncogenicExpected
                    
        listOfDicts.append({'Hugo_Symbol': gene, 'N_observed': curNumObserved, 'N_expected': curNumExpected})
    return pd.DataFrame(listOfDicts)



def summarize_observed_and_expected_oncogenic_mut_rates(maf, oncogenicSDict, sigs, genes=set(['TP53'])):
    oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
    impactVersionDict = {'IM3': 'IMPACT_368', 'IM5': 'IMPACT_410', 'IM6': 'IMPACT_468'}
    panelGeneDict = {
            'IMPACT_368': set(['ABL1', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BAP1', 'BARD1', 'BBC3', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CARD11', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79B', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSF1R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EP300', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERG', 'ESR1', 'ETV1', 'ETV6', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXP1', 'FUBP1', 'GATA1', 'GATA2', 'GATA3', 'GNA11', 'GNAQ', 'GNAS', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3B', 'HNF1A', 'HRAS', 'ICOSLG', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INPP4A', 'INPP4B', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAPK1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH6', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOR1', 'NF1', 'NF2', 'NFE2L2', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTRK1', 'NTRK2', 'NTRK3', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLE', 'PPP2R1A', 'PRDM1', 'PRKAR1A', 'PTCH1', 'PTEN', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAC1', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RUNX1', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SETD2', 'SF3B1', 'SH2D1A', 'SHQ1', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SOCS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SRC', 'STAG2', 'STK11', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TBX3', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP63', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'VHL', 'VTCN1', 'WT1', 'XIAP', 'XPO1', 'YAP1', 'YES1']),
            'IMPACT_410': set(['ABL1', 'ACVR1', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'EP300', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH6', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPM1D', 'PPP2R1A', 'PPP6C', 'PRDM1', 'PRKAR1A', 'PTCH1', 'PTEN', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAD21', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RUNX1', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHQ1', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SOCS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'VEGFA', 'VHL', 'VTCN1', 'WT1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2']),
            'IMPACT_468': set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'AMER1', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'GTF2I', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MLL4', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SETD8', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2']) 
            }
    cntr = 0
    listOfDicts = []
    for case in set(maf['Tumor_Sample_Barcode']):
        versionKey = case[-3:]
        impactVersion = None
        if versionKey in impactVersionDict:
            impactVersion = impactVersionDict[versionKey]
        
        cntr += 1
        if cntr%10 == 0: print cntr
        
        if impactVersion != None:
            caseSigs = sigs[sigs['Tumor_Sample_Barcode'] == case]
            caseMaf = maf[maf['Tumor_Sample_Barcode'] == case]
            genes = panelGeneDict[impactVersion]
            geneMuts = caseMaf[caseMaf['Hugo_Symbol'].isin(genes)]
            oncogenicMuts = geneMuts[geneMuts['oncogenic'].isin(oncogenicMutColNames)]
            nOncogenicSnps = oncogenicMuts[oncogenicMuts['Variant_Type'] == 'SNP'].shape[0]
            nSnps = geneMuts[geneMuts['Variant_Type'] == 'SNP'].shape[0]
            nMuts = geneMuts.shape[0]
            if caseSigs.shape[0] > 0:
                signatureColumnNames = ['mean_' + str(i) for i in range(1,31)]
                signatureColsOnly = caseSigs[signatureColumnNames]
                renameDict = dict([('mean_' + str(i), 'Signature.' + str(i)) for i in range(1,31)])
                signatureColsOnly = signatureColsOnly.rename(columns=renameDict)
                
                decompositionDict = signatureColsOnly.to_dict(orient='records')[0]
                quadNucFractions = mutation_modeling_util.get_quadnuc_fracs_given_decomposition(decompositionDict, spectraPath = pathPrefix + '/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt')
                oncogenicExpected = mutation_modeling_util.get_expected_oncogenic_val_given_quadnuc_fractions(quadNucFractions, oncogenicSDict, impactVersion)
                nOncogenicExpected = oncogenicExpected * nSnps
            
                
                listOfDicts.append({'Tumor_Sample_Barcode': case, 'Nmut': nMuts, 'NSnps': nSnps, 'NOncogenicSnpsObserved': nOncogenicSnps, 'nOncogenicSnpsExpected': nOncogenicExpected})
    
    return pd.DataFrame(listOfDicts)


def get_limited_observed_vs_expected_from_cohort(mMaf, ids, gSpecificOncSusceptibilityDict, sigs, domSig, n=15, t= .1):
    caseMaf = mMaf[mMaf['Tumor_Sample_Barcode'].isin(ids)]
    l = get_top_n_most_mutable_genes_by_signature(domSig, geneSpecificOncogenicSusceptibilityDict, n)
    mutableNames = set([i[0] for i in l])
    frequentlyMutatedNames = get_genes_mutated_in_cohort_above_threshold(caseMaf, thresh=t)
    genesToAnalyze = mutableNames | frequentlyMutatedNames
    print genesToAnalyze
    observedAndExpectedGeneInfoDf = compare_observed_vs_expected_gene_oncogenic_mut_prevalence_across_cohort(caseMaf, gSpecificOncSusceptibilityDict, sigs, genes=genesToAnalyze)
    return observedAndExpectedGeneInfoDf 
 

#TODO ALERT!!!!! MOVE THIS TO ANOTHER SCRIPT AREA
def summarize_cohort_by_gene(mMaf, hyperIds, normalIds):
    
    d = {}
    
    oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
    caseMafHyper = mMaf[mMaf['Tumor_Sample_Barcode'].isin(hyperIds)]
    caseMafNormal = mMaf[mMaf['Tumor_Sample_Barcode'].isin(normalIds)]
    hyperOncogenic = caseMafHyper[caseMafHyper['oncogenic'].isin(oncogenicMutColNames)]
    normalOncogenic = caseMafNormal[caseMafNormal['oncogenic'].isin(oncogenicMutColNames)]
    nHyper = len(hyperIds)
    nNormal = len(normalIds)
    
    cntr = 0
    for gene in set(caseMafHyper['Hugo_Symbol']):
        if cntr%10 == 0:
            print cntr
        cntr += 1
        geneMafNormal = normalOncogenic[normalOncogenic['Hugo_Symbol'] == gene]
        geneMafHyper = hyperOncogenic[hyperOncogenic['Hugo_Symbol'] == gene]
        fracNormal = 1.0*geneMafNormal.shape[0]/nNormal
        fracHyperHotspot = 0
        if geneMafHyper.shape[0] != 0:
            fracHyperHotspot = 1.0*geneMafHyper[geneMafHyper['is-a-hotspot'] == 'Y'].shape[0]/geneMafHyper.shape[0]
        d[gene] = (fracNormal, fracHyperHotspot)
    return d
        
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

def enumerate_observed_expected_info_for_cancer_type(mutationMaf, cancerType, signature, geneSpecificOncogenicSusceptibilityDict, impactSigs, cancerTypeSimData=None):
    cancerTypeHyperIds = set(impactSigs[(impactSigs['Nmut_Mb'] > 50) & (impactSigs['cancer_type'] == cancerType)]['Tumor_Sample_Barcode'])
    cancerTypeNotHyperIds = set(impactSigs[(impactSigs['Nmut_Mb'] < 50) & (impactSigs['cancer_type'] == cancerType)]['Tumor_Sample_Barcode'])
    if cancerTypeSimData == None:
        cancerTypeSimData = get_limited_observed_vs_expected_from_cohort(mutationMaf, cancerTypeHyperIds, geneSpecificOncogenicSusceptibilityDict, impactSigs, signature)
    cancerTypeIncidenceDict = summarize_cohort_by_gene(mutationMaf, cancerTypeHyperIds, cancerTypeNotHyperIds)

    cancerTypeSimData = cancerTypeSimData[cancerTypeSimData['N_expected'] > 0]
    nCases = len(curIds)
    tumorSupressors = set(['ERRFI1', 'ASXL2', 'PMAIP1', 'ACTG1', 'SUFU', 'FBXO11', 'MEN1', 'FAM58A', 'B2M', 'RB1', 'DUSP22', 'SESN1', 'GPS2', 'RAD51D', 'SMG1', 'CDC73', 'MAP3K1', 'SMARCB1', 'INPP4B', 'PARK2', 'SMAD4', 'CBFB', 'CDH1', 'PPP6C', 'SETDB1', 'SETDB2', 'NF2', 'CDKN2B', 'CDKN2C', 'CDKN2A', 'DDX3X', 'PIK3R1', 'BARD1', 'PDS5B', 'KLF4', 'SPRED1', 'VHL', 'SMAD2', 'PMS1', 'PMS2', 'SETD2', 'GATA3', 'TBL1XR1', 'MUTYH', 'SOCS1', 'FAM175A', 'ROBO1', 'ARID1B', 'ARID1A', 'TCF7L2', 'STK11', 'FOXA1', 'PTEN', 'FAT1', 'FAS', 'CYLD', 'MAX', 'SH2D1A', 'APC', 'NTHL1', 'CTCF', 'KDM5C', 'KMT2C', 'ZFHX3', 'FOXP1', 'PIGA', 'CDKN1B', 'CDKN1A', 'FUBP1', 'MSH2', 'ID3', 'TNFRSF14', 'TRAF3', 'EP400', 'BRIP1', 'ARID4A', 'ARID4B', 'XRCC2', 'DAXX', 'SDHAF2', 'ASXL1', 'AMER1', 'RASA1', 'EGR1', 'MST1', 'SOX17', 'RUNX1', 'PIK3R3', 'NCOR1', 'NF1', 'JAK1', 'PTPRD', 'CHEK2', 'CHEK1', 'SMC1A', 'TMEM127', 'STAG1', 'RAD51', 'TCF3', 'STAG2', 'ARID2', 'RAD50', 'RNF43', 'PARP1', 'BLM', 'CUX1', 'RECQL', 'RAD21', 'PTPN2', 'PTPN1', 'SLX4', 'INHA', 'PAX5', 'IRF1', 'TP53', 'HLA-A', 'IRF8', 'CBL', 'TOP1', 'SHQ1', 'PRDM1', 'NSD1', 'ATXN2', 'CREBBP', 'HDAC4', 'SESN2', 'PPP2R1A', 'EPHA7', 'ATM', 'EPHA3', 'POT1', 'SMAD3', 'MOB3B', 'TBX3', 'POLE', 'ATR', 'FANCD2', 'FH', 'BCORL1', 'SOX9', 'IKZF3', 'TSC1', 'TP63', 'MRE11A', 'SDHC', 'BTG1', 'POLD1', 'CIITA', 'SMC3', 'SAMHD1', 'RTEL1', 'ECT2L', 'PIK3R2', 'CRBN', 'FANCC', 'NBN', 'FANCA', 'HLA-B', 'RECQL4', 'DUSP4', 'ERCC2', 'FBXW7', 'TGFBR2', 'TGFBR1', 'MSH3', 'RBM15', 'TET1', 'TET3', 'SESN3', 'MGA', 'LTB', 'FOXL2', 'SH2B3', 'BCOR', 'HIST1H1D', 'ATRX', 'EP300', 'RAD51C', 'RAD51B', 'HIST1H1B', 'TNFAIP3', 'DICER1', 'ARID5B', 'LATS2', 'FOXO1', 'KEAP1', 'EZH2', 'SP140', 'NKX3-1', 'PBRM1', 'PALB2', 'CIC', 'BRCA1', 'DTX1', 'FLCN', 'SPEN', 'CD58', 'ERCC3', 'ERCC4', 'MSH6', 'BCL11B', 'BMPR1A', 'ERF', 'BRCA2', 'NOTCH2', 'EED', 'MITF', 'ELF3', 'SMARCA4', 'BBC3', 'ANKRD11', 'CEBPA', 'BCL2L11', 'AXIN2', 'AXIN1', 'CDK12', 'ESCO2', 'MLH1', 'SDHB', 'MED12', 'HNF1A', 'RYBP', 'ATP6V1B2', 'DNMT3B', 'KMT2B', 'KMT2A', 'DNMT3A', 'NFKBIA', 'TRAF5', 'KMT2D', 'SPOP', 'RBM10', 'P2RY8', 'TP53BP1', 'TSC2', 'KDM6A', 'EPCAM', 'PHOX2B', 'NPM1', 'BCL10', 'LATS1', 'HOXB13', 'ARID3A', 'PTPRT', 'PTPRS', 'INPPL1', 'NOTCH4', 'TET2', 'NOTCH1', 'CASP8', 'NOTCH3', 'GRIN2A', 'MAP2K4', 'WT1', 'BACH2', 'SDHA', 'BAP1', 'PTCH1', 'SDHD'])
    cancerTypeSimData = cancerTypeSimData[cancerTypeSimData['Hugo_Symbol'].isin(tumorSupressors)]
    
    recurrentTumorSupressors, recurrentOncogenes = enumerate_recurrently_mutated_tumor_supressors_and_oncogenes(mutationMaf[mutationMaf['Tumor_Sample_Barcode'].isin(cancerTypeNotHyperIds)])
    
    cancerTypeSimData['observedMutPerCase'] = cancerTypeSimData['N_observed']/nCases
    cancerTypeSimData['expectedMutPerCase'] = cancerTypeSimData['N_expected']/nCases
    cancerTypeSimData['normalCohortFrac'] = cancerTypeSimData['Hugo_Symbol'].apply(lambda x: cancerTypeIncidenceDict[x][0] if x in cancerTypeIncidenceDict else None)
    cancerTypeSimData['hyperHotspotFrac'] = cancerTypeSimData['Hugo_Symbol'].apply(lambda x: cancerTypeIncidenceDict[x][1] if x in cancerTypeIncidenceDict else None)
    
    writeFileName = '~/Desktop/WORK/dataForLocalPlotting/geneSpecificSimulationData_' + cancerType + '.tsv'
    print 'writing to ', writeFileName
    cancerTypeSimData.to_csv(writeFileName, sep='\t', index=False)
    return cancerTypeSimData

     
"""validQuadNucs = set([start + change + end for start in ['A', 'C', 'T', 'G'] 
                        for change in ['CA', 'CG', 'CT', 'TA', 'TC', 'TG']
                        for end in ['A', 'C', 'T', 'G']])"""

#CODE to look at all possible impact muts and combine stuff together
#mutData = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/simulatedMafs/allmut_impact_simulated_muts_mafAnno_trinuc.maf') 
geneDataSummaryDf = mutation_modeling_util.load_in_and_summarize_simulated_mut_information(mafDirPath = pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/simulatedMafs')
geneDataSummaryDf.to_csv(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/panImpactSimulatedMutationDataSummary.tsv', sep='\t', index=False)

possibleMutData = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/panImpactSimulatedMutationDataSummary.tsv')

mutationMaf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/hypermutationAnalysisProjectData/hypermutationAnalysisProjectMaf.tsv')
mutationMaf['quadNuc'] = mutationMaf.apply(lambda row: mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)

impactSigs = pd.read_table(pathPrefix + '/ifs/res/taylorlab/impact_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
impactSigs['pid'] = impactSigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
impactSigs['cancer_type'] = impactSigs['pid'].apply(lambda x: cDict[x] if x in cDict else None)

impactSigs.to_csv('~/Desktop/WORK/dataForLocalPlotting/vanillaImpactSigsFile.tsv', sep='\t', index=False)










#TODO adjust the possible mut data so we can reason on it specific to each panel
oncogenicSusceptibilityDict = mutation_modeling_util.calculate_quadnuc_based_oncogenic_susceptibility_dict(possibleMutData)
geneSpecificOncogenicSusceptibilityDict = mutation_modeling_util.calculate_quadnuc_and_gene_based_oncogenic_susceptibility_dict(possibleMutData)
#LETS start by calculating the differential oncogenic rates by signature
df = calculate_oncogenic_mut_susceptibility_of_genes_by_signature(oncogenicSusceptibilityDict)
df.to_csv('~/Desktop/WORK/dataForLocalPlotting/oncogenicSusceptibilityBySignature.tsv', sep='\t', index=False)

##((((((((#####################(################)))))))))
#NEXT we are doing the expected vs observed mut rates
dataF = summarize_observed_and_expected_oncogenic_mut_rates(mutationMaf, oncogenicSusceptibilityDict, impactSigs)
#dataF['ratio'] = dataF.apply(lambda row: None if row['NOncogenicObserved'] == 0 else row['NOncogenicObserved']/row['nOncogenicExpected'], axis=1)
#dataF['ratio'] = dataF.apply(lambda row: None if row['NOncogenicObserved'] == 0 else math.log(row['NOncogenicObserved']/row['nOncogenicSnpsExpected']), axis=1)

#TODO change this based on SNP burden TODO TODO TODO TODO 
dataF['percentile'] = pd.qcut(dataF['NSnps'], [0, .5, .75, .85, .9, .95, .97, .99, 1], 
     labels=['0-50th', '50-75th', '75th-85th', '85th-90th', '90th-95th', '95th-97th', '97th-99th', '99th'])

dataF['percentile2'] = pd.qcut(dataF['Nmut'], 
     [0,.25,.5,.6,.7, .8] + [.8 + .01*i for i in range(1,21)],
     ['0-25th', '25-50th', '50-60th', '60-70th', '70th-80th', '80th'] + [str(80 + i) + 'th' for i in range(1,20)])

dataF.to_csv('~/Desktop/WORK/dataForLocalPlotting/simulationData.tsv', index=False, sep='\t')

labs = ['0-50th', '50-75th', '75th-85th', '85th-90th', '90th-95th', '95th-97th', '97th-99th', '99th']
for lab in labs:
    print lab, np.nanmedian(dataF[dataF['percentile'] == lab]['NOncogenicSnpsObserved']) - np.nanmedian(dataF[dataF['percentile'] == lab]['nOncogenicSnpsExpected'])

#%#%#%#%#%#%#%%#
#NEXT LOOK AT PER GENE MUTATION OBSERVATIONS

#l = get_top_n_most_mutable_genes_by_signature('Signature.1', geneSpecificOncogenicSusceptibilityDict, 15)

tumorSupressors = set(['ERRFI1', 'ASXL2', 'PMAIP1', 'ACTG1', 'SUFU', 'FBXO11', 'MEN1', 'FAM58A', 'B2M', 'RB1', 'DUSP22', 'SESN1', 'GPS2', 'RAD51D', 'SMG1', 'CDC73', 'MAP3K1', 'SMARCB1', 'INPP4B', 'PARK2', 'SMAD4', 'CBFB', 'CDH1', 'PPP6C', 'SETDB1', 'SETDB2', 'NF2', 'CDKN2B', 'CDKN2C', 'CDKN2A', 'DDX3X', 'PIK3R1', 'BARD1', 'PDS5B', 'KLF4', 'SPRED1', 'VHL', 'SMAD2', 'PMS1', 'PMS2', 'SETD2', 'GATA3', 'TBL1XR1', 'MUTYH', 'SOCS1', 'FAM175A', 'ROBO1', 'ARID1B', 'ARID1A', 'TCF7L2', 'STK11', 'FOXA1', 'PTEN', 'FAT1', 'FAS', 'CYLD', 'MAX', 'SH2D1A', 'APC', 'NTHL1', 'CTCF', 'KDM5C', 'KMT2C', 'ZFHX3', 'FOXP1', 'PIGA', 'CDKN1B', 'CDKN1A', 'FUBP1', 'MSH2', 'ID3', 'TNFRSF14', 'TRAF3', 'EP400', 'BRIP1', 'ARID4A', 'ARID4B', 'XRCC2', 'DAXX', 'SDHAF2', 'ASXL1', 'AMER1', 'RASA1', 'EGR1', 'MST1', 'SOX17', 'RUNX1', 'PIK3R3', 'NCOR1', 'NF1', 'JAK1', 'PTPRD', 'CHEK2', 'CHEK1', 'SMC1A', 'TMEM127', 'STAG1', 'RAD51', 'TCF3', 'STAG2', 'ARID2', 'RAD50', 'RNF43', 'PARP1', 'BLM', 'CUX1', 'RECQL', 'RAD21', 'PTPN2', 'PTPN1', 'SLX4', 'INHA', 'PAX5', 'IRF1', 'TP53', 'HLA-A', 'IRF8', 'CBL', 'TOP1', 'SHQ1', 'PRDM1', 'NSD1', 'ATXN2', 'CREBBP', 'HDAC4', 'SESN2', 'PPP2R1A', 'EPHA7', 'ATM', 'EPHA3', 'POT1', 'SMAD3', 'MOB3B', 'TBX3', 'POLE', 'ATR', 'FANCD2', 'FH', 'BCORL1', 'SOX9', 'IKZF3', 'TSC1', 'TP63', 'MRE11A', 'SDHC', 'BTG1', 'POLD1', 'CIITA', 'SMC3', 'SAMHD1', 'RTEL1', 'ECT2L', 'PIK3R2', 'CRBN', 'FANCC', 'NBN', 'FANCA', 'HLA-B', 'RECQL4', 'DUSP4', 'ERCC2', 'FBXW7', 'TGFBR2', 'TGFBR1', 'MSH3', 'RBM15', 'TET1', 'TET3', 'SESN3', 'MGA', 'LTB', 'FOXL2', 'SH2B3', 'BCOR', 'HIST1H1D', 'ATRX', 'EP300', 'RAD51C', 'RAD51B', 'HIST1H1B', 'TNFAIP3', 'DICER1', 'ARID5B', 'LATS2', 'FOXO1', 'KEAP1', 'EZH2', 'SP140', 'NKX3-1', 'PBRM1', 'PALB2', 'CIC', 'BRCA1', 'DTX1', 'FLCN', 'SPEN', 'CD58', 'ERCC3', 'ERCC4', 'MSH6', 'BCL11B', 'BMPR1A', 'ERF', 'BRCA2', 'NOTCH2', 'EED', 'MITF', 'ELF3', 'SMARCA4', 'BBC3', 'ANKRD11', 'CEBPA', 'BCL2L11', 'AXIN2', 'AXIN1', 'CDK12', 'ESCO2', 'MLH1', 'SDHB', 'MED12', 'HNF1A', 'RYBP', 'ATP6V1B2', 'DNMT3B', 'KMT2B', 'KMT2A', 'DNMT3A', 'NFKBIA', 'TRAF5', 'KMT2D', 'SPOP', 'RBM10', 'P2RY8', 'TP53BP1', 'TSC2', 'KDM6A', 'EPCAM', 'PHOX2B', 'NPM1', 'BCL10', 'LATS1', 'HOXB13', 'ARID3A', 'PTPRT', 'PTPRS', 'INPPL1', 'NOTCH4', 'TET2', 'NOTCH1', 'CASP8', 'NOTCH3', 'GRIN2A', 'MAP2K4', 'WT1', 'BACH2', 'SDHA', 'BAP1', 'PTCH1', 'SDHD'])

endoSimData = enumerate_observed_expected_info_for_cancer_type(mutationMaf, 'Endometrial Cancer', 'Signature.10', geneSpecificOncogenicSusceptibilityDict, impactSigs)
endoSimDataLim = endoSimData[endoSimData['Hugo_Symbol'].isin(tumorSupressors)]
recurrentTumorSupressors, recurrentOncogenes = enumerate_recurrently_mutated_tumor_supressors_and_oncogenes(mutationMaf[mutationMaf['Tumor_Sample_Barcode'].isin(endometrialNotHyperIds)])
endoSimDataLim['DisplayName'] =  endoSimDataLim.apply(lambda row: row['Hugo_Symbol'] if row['observedMutPerCase'] > .25 else None, axis=1)
endoSimDataLim['isRecurrent'] =  endoSimDataLim['Hugo_Symbol'].apply(lambda x: True if x in recurrentTumorSupressors else False)
writeFileName = '~/Desktop/WORK/dataForLocalPlotting/geneSpecificSimulationData_' + 'Endometrial Cancer' + '.tsv'
endoSimDataLim.to_csv(writeFileName, sep='\t', index=False)

colorectalSimData = enumerate_observed_expected_info_for_cancer_type(mutationMaf, 'Colorectal Cancer', 'Signature.6', geneSpecificOncogenicSusceptibilityDict, impactSigs)
colorectalSimDataLim = colorectalSimData[colorectalSimData['Hugo_Symbol'].isin(tumorSupressors)]
recurrentTumorSupressors, recurrentOncogenes = enumerate_recurrently_mutated_tumor_supressors_and_oncogenes(mutationMaf[mutationMaf['Tumor_Sample_Barcode'].isin(colorectalNotHyperIds)])
colorectalSimDataLim['DisplayName'] =  colorectalSimDataLim.apply(lambda row: row['Hugo_Symbol'] if row['observedMutPerCase'] > .25 else row['Hugo_Symbol'] if row['Hugo_Symbol'] == 'RNF43' or row['expectedMutPerCase'] > .3 else None, axis=1)
colorectalSimDataLim['isRecurrent'] =  colorectalSimDataLim['Hugo_Symbol'].apply(lambda x: True if x in recurrentTumorSupressors else False)
writeFileName = '~/Desktop/WORK/dataForLocalPlotting/geneSpecificSimulationData_' + 'Colorectal_Cancer' + '.tsv'
colorectalSimDataLim.to_csv(writeFileName, sep='\t', index=False)


gliomaSimData = enumerate_observed_expected_info_for_cancer_type(mutationMaf, 'Glioma', 'Signature.11', geneSpecificOncogenicSusceptibilityDict, impactSigs)
gliomaSimDataLim = gliomaSimData[gliomaSimData['Hugo_Symbol'].isin(tumorSupressors)]
recurrentTumorSupressors, recurrentOncogenes = enumerate_recurrently_mutated_tumor_supressors_and_oncogenes(mutationMaf[mutationMaf['Tumor_Sample_Barcode'].isin(gliomaNotHyperIds)])
gliomaSimDataLim['DisplayName'] =  gliomaSimDataLim.apply(lambda row: row['Hugo_Symbol'] if row['observedMutPerCase'] > .2 else row['Hugo_Symbol'] if row['expectedMutPerCase'] > .2 
                else row['Hugo_Symbol'] if row['Hugo_Symbol'] == 'MLH1' or row['Hugo_Symbol'] == 'MSH2' or row['Hugo_Symbol'] == 'PMS2' else  None, axis=1)
gliomaSimDataLim['isRecurrent'] =  gliomaSimDataLim['Hugo_Symbol'].apply(lambda x: True if x in recurrentTumorSupressors else False)
writeFileName = '~/Desktop/WORK/dataForLocalPlotting/geneSpecificSimulationData_' + 'glioma' + '.tsv'
gliomaSimDataLim.to_csv(writeFileName, sep='\t', index=False)


gliomaNotHyperIds = set(impactSigs[(impactSigs['Nmut_Mb'] < 50) & (impactSigs['cancer_type'] == 'Glioma')]['Tumor_Sample_Barcode'])


cdsLengthData = pd.read_table('~/Desktop')


"""endometrialHyperIds = set(impactSigs[(impactSigs['Nmut_Mb'] > 50) & (impactSigs['cancer_type'] == 'Endometrial Cancer')]['Tumor_Sample_Barcode'])
endometrialNotHyperIds = set(impactSigs[(impactSigs['Nmut_Mb'] < 50) & (impactSigs['cancer_type'] == 'Endometrial Cancer')]['Tumor_Sample_Barcode'])
endometrialSimData = get_limited_observed_vs_expected_from_cohort(mutationMaf, endometrialHyperIds, geneSpecificOncogenicSusceptibilityDict, impactSigs, 'Signature.10')

colorectalHyperIds = set(impactSigs[(impactSigs['Nmut_Mb'] > 50) & (impactSigs['cancer_type'] == 'Colorectal Cancer')]['Tumor_Sample_Barcode'])
colorectalNotHyperIds = set(impactSigs[(impactSigs['Nmut_Mb'] < 50) & (impactSigs['cancer_type'] == 'Colorectal Cancer')]['Tumor_Sample_Barcode'])
colorectalSimData = get_limited_observed_vs_expected_from_cohort(mutationMaf, colorectalHyperIds, geneSpecificOncogenicSusceptibilityDict, impactSigs, 'Signature.6')

gliomaHyperIds = set(impactSigs[(impactSigs['Nmut_Mb'] > 50) & (impactSigs['cancer_type'] == 'Glioma')]['Tumor_Sample_Barcode'])
gliomaNotHyperIds = set(impactSigs[(impactSigs['Nmut_Mb'] < 50) & (impactSigs['cancer_type'] == 'Glioma')]['Tumor_Sample_Barcode'])
gliomaSimData = get_limited_observed_vs_expected_from_cohort(mutationMaf, gliomaHyperIds, geneSpecificOncogenicSusceptibilityDict, impactSigs, 'Signature.11')

bladderHyperIds = set(impactSigs[(impactSigs['Nmut_Mb'] > 50) & (impactSigs['cancer_type'] == 'Bladder Cancer')]['Tumor_Sample_Barcode'])
bladderNotHyperIds = set(impactSigs[(impactSigs['Nmut_Mb'] < 50) & (impactSigs['cancer_type'] == 'Bladder Cancer')]['Tumor_Sample_Barcode'])
bladderSimData = get_limited_observed_vs_expected_from_cohort(mutationMaf, bladderHyperIds, geneSpecificOncogenicSusceptibilityDict, impactSigs, 'Signature.2')

melanomaHyperIds = set(impactSigs[(impactSigs['Nmut_Mb'] > 50) & (impactSigs['cancer_type'] == 'Melanoma')]['Tumor_Sample_Barcode'])
melanomaNotHyperIds = set(impactSigs[(impactSigs['Nmut_Mb'] < 50) & (impactSigs['cancer_type'] == 'Melanoma')]['Tumor_Sample_Barcode'])
melanomaSimData = get_limited_observed_vs_expected_from_cohort(mutationMaf, melanomaHyperIds, geneSpecificOncogenicSusceptibilityDict, impactSigs, 'Signature.7')


endometrialIncidenceDict = summarize_cohort_by_gene(mutationMaf, endometrialHyperIds, endometrialNotHyperIds)
colorectalIncidenceDict = summarize_cohort_by_gene(mutationMaf, colorectalHyperIds, colorectalNotHyperIds)
gliomaIncidenceDict = summarize_cohort_by_gene(mutationMaf, gliomaHyperIds, gliomaNotHyperIds)
bladderIncidenceDict = summarize_cohort_by_gene(mutationMaf, bladderHyperIds, bladderNotHyperIds)
melanomaIncidenceDict = summarize_cohort_by_gene(mutationMaf, melanomaHyperIds, melanomaNotHyperIds)


curIds = endometrialHyperIds
incidenceDict = endometrialIncidenceDict
caseSimData = endometrialSimData

caseSimData = caseSimData[caseSimData['N_expected'] > 0]
nCases = len(curIds)
tumorSupressors = set(['ERRFI1', 'ASXL2', 'PMAIP1', 'ACTG1', 'SUFU', 'FBXO11', 'MEN1', 'FAM58A', 'B2M', 'RB1', 'DUSP22', 'SESN1', 'GPS2', 'RAD51D', 'SMG1', 'CDC73', 'MAP3K1', 'SMARCB1', 'INPP4B', 'PARK2', 'SMAD4', 'CBFB', 'CDH1', 'PPP6C', 'SETDB1', 'SETDB2', 'NF2', 'CDKN2B', 'CDKN2C', 'CDKN2A', 'DDX3X', 'PIK3R1', 'BARD1', 'PDS5B', 'KLF4', 'SPRED1', 'VHL', 'SMAD2', 'PMS1', 'PMS2', 'SETD2', 'GATA3', 'TBL1XR1', 'MUTYH', 'SOCS1', 'FAM175A', 'ROBO1', 'ARID1B', 'ARID1A', 'TCF7L2', 'STK11', 'FOXA1', 'PTEN', 'FAT1', 'FAS', 'CYLD', 'MAX', 'SH2D1A', 'APC', 'NTHL1', 'CTCF', 'KDM5C', 'KMT2C', 'ZFHX3', 'FOXP1', 'PIGA', 'CDKN1B', 'CDKN1A', 'FUBP1', 'MSH2', 'ID3', 'TNFRSF14', 'TRAF3', 'EP400', 'BRIP1', 'ARID4A', 'ARID4B', 'XRCC2', 'DAXX', 'SDHAF2', 'ASXL1', 'AMER1', 'RASA1', 'EGR1', 'MST1', 'SOX17', 'RUNX1', 'PIK3R3', 'NCOR1', 'NF1', 'JAK1', 'PTPRD', 'CHEK2', 'CHEK1', 'SMC1A', 'TMEM127', 'STAG1', 'RAD51', 'TCF3', 'STAG2', 'ARID2', 'RAD50', 'RNF43', 'PARP1', 'BLM', 'CUX1', 'RECQL', 'RAD21', 'PTPN2', 'PTPN1', 'SLX4', 'INHA', 'PAX5', 'IRF1', 'TP53', 'HLA-A', 'IRF8', 'CBL', 'TOP1', 'SHQ1', 'PRDM1', 'NSD1', 'ATXN2', 'CREBBP', 'HDAC4', 'SESN2', 'PPP2R1A', 'EPHA7', 'ATM', 'EPHA3', 'POT1', 'SMAD3', 'MOB3B', 'TBX3', 'POLE', 'ATR', 'FANCD2', 'FH', 'BCORL1', 'SOX9', 'IKZF3', 'TSC1', 'TP63', 'MRE11A', 'SDHC', 'BTG1', 'POLD1', 'CIITA', 'SMC3', 'SAMHD1', 'RTEL1', 'ECT2L', 'PIK3R2', 'CRBN', 'FANCC', 'NBN', 'FANCA', 'HLA-B', 'RECQL4', 'DUSP4', 'ERCC2', 'FBXW7', 'TGFBR2', 'TGFBR1', 'MSH3', 'RBM15', 'TET1', 'TET3', 'SESN3', 'MGA', 'LTB', 'FOXL2', 'SH2B3', 'BCOR', 'HIST1H1D', 'ATRX', 'EP300', 'RAD51C', 'RAD51B', 'HIST1H1B', 'TNFAIP3', 'DICER1', 'ARID5B', 'LATS2', 'FOXO1', 'KEAP1', 'EZH2', 'SP140', 'NKX3-1', 'PBRM1', 'PALB2', 'CIC', 'BRCA1', 'DTX1', 'FLCN', 'SPEN', 'CD58', 'ERCC3', 'ERCC4', 'MSH6', 'BCL11B', 'BMPR1A', 'ERF', 'BRCA2', 'NOTCH2', 'EED', 'MITF', 'ELF3', 'SMARCA4', 'BBC3', 'ANKRD11', 'CEBPA', 'BCL2L11', 'AXIN2', 'AXIN1', 'CDK12', 'ESCO2', 'MLH1', 'SDHB', 'MED12', 'HNF1A', 'RYBP', 'ATP6V1B2', 'DNMT3B', 'KMT2B', 'KMT2A', 'DNMT3A', 'NFKBIA', 'TRAF5', 'KMT2D', 'SPOP', 'RBM10', 'P2RY8', 'TP53BP1', 'TSC2', 'KDM6A', 'EPCAM', 'PHOX2B', 'NPM1', 'BCL10', 'LATS1', 'HOXB13', 'ARID3A', 'PTPRT', 'PTPRS', 'INPPL1', 'NOTCH4', 'TET2', 'NOTCH1', 'CASP8', 'NOTCH3', 'GRIN2A', 'MAP2K4', 'WT1', 'BACH2', 'SDHA', 'BAP1', 'PTCH1', 'SDHD'])
caseSimData = caseSimData[caseSimData['Hugo_Symbol'].isin(tumorSupressors)]
    
recurrentTumorSupressors, recurrentOncogenes = enumerate_recurrently_mutated_tumor_supressors_and_oncogenes(mutationMaf[mutationMaf['Tumor_Sample_Barcode'].isin(endometrialNotHyperIds)])

#geneNames = set([''])
#caseSimData['displayName'] = caseSimData.apply(lambda row: row['Hugo_Symbol'] if , axis=1)
caseSimData['observedMutPerCase'] = caseSimData['N_observed']/nCases
caseSimData['expectedMutPerCase'] = caseSimData['N_expected']/nCases
caseSimData['normalCohortFrac'] = caseSimData['Hugo_Symbol'].apply(lambda x: incidenceDict[x][0] if x in incidenceDict else None)
caseSimData['hyperHotspotFrac'] = caseSimData['Hugo_Symbol'].apply(lambda x: incidenceDict[x][1] if x in incidenceDict else None)
caseSimData.to_csv('~/Desktop/WORK/dataForLocalPlotting/geneSpecificSimulationData_endometrial.tsv', sep='\t', index=False)
"""










