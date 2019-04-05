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
import mutation_modeling_util


    
    
    
def model_oncogenic_mut_acquisition(allSpectraDict, qNucDict, qNucDictAdj, geneOncogenicOrHotspotDict, geneMutProbInfo, simulatedNmuts = 1000, simulatedKCases = 20, mutationMode='oncogenic'):
    
    #initiate the signatures dicts
    signatureCols = ['Signature.' + str(i) for i in range(1,31)]
    sigsDict = dict()
    for col in signatureCols:
        sigsDict[col] = 0
     
    simulatedNmuts = 1000
    simulatedKCases = 20
    simulatedVals = []
    for col in signatureCols:
        sigsDict[col] = 1
        val1 = mutation_modeling_util.do_k_simulations_of_mutations(simulatedKCases,                                      
                                      sigsDict, 
                                      allSpectraDict, #a dictionary of all spectras from which we get a mutation signature spectra
                                      simulatedNmuts, #the number of mutations we are supposed to simulate
                                      qNucDict, #probability of picking a specific quadnuc given a spectra 
                                      qNucDictAdj,
                                      geneOncogenicOrHotspotDict, #probability of an oncogenic mutation O given a mutation in a gene G and a quadnuc Q
                                      geneMutProbInfo, #The probability of mutating a gene (based on the probabilities of genes being mutated across all silent mutations)
                                      mutMode=mutationMode, #mutMode--do we care about picking oncogenic mutations or picking hotspot mutations
                                      modelName='model1' #model: specifies the assumptions we use to pick the mutations
                                      )
        
        val2 = mutation_modeling_util.do_k_simulations_of_mutations(simulatedKCases,                                      
                                      sigsDict, 
                                      allSpectraDict, #a dictionary of all spectras from which we get a mutation signature spectra
                                      simulatedNmuts, #the number of mutations we are supposed to simulate
                                      qNucDict, #probability of picking a specific quadnuc given a spectra 
                                      qNucDictAdj,
                                      geneOncogenicOrHotspotDict, #probability of an oncogenic mutation O given a mutation in a gene G and a quadnuc Q
                                      geneMutProbInfo, #The probability of mutating a gene (based on the probabilities of genes being mutated across all silent mutations)
                                      mutMode=mutationMode, #mutMode--do we care about picking oncogenic mutations or picking hotspot mutations
                                      modelName='model2' #model: specifies the assumptions we use to pick the mutations
                                      )
        
        val3 = mutation_modeling_util.do_k_simulations_of_mutations(simulatedKCases,                                      
                                      sigsDict, 
                                      allSpectraDict, #a dictionary of all spectras from which we get a mutation signature spectra
                                      simulatedNmuts, #the number of mutations we are supposed to simulate
                                      qNucDict, #probability of picking a specific quadnuc given a spectra 
                                      qNucDictAdj,
                                      geneOncogenicOrHotspotDict, #probability of an oncogenic mutation O given a mutation in a gene G and a quadnuc Q
                                      geneMutProbInfo, #The probability of mutating a gene (based on the probabilities of genes being mutated across all silent mutations)
                                      mutMode=mutationMode, #mutMode--do we care about picking oncogenic mutations or picking hotspot mutations
                                      modelName='model3' #model: specifies the assumptions we use to pick the mutations
                                      )
        
        simulatedVals.append((col, val1, val2, val3))
        sigsDict[col] = 0
        
    return simulatedVals

def convert_dict_to_df(d):
    genes = d.keys()
    listOfDicts = []
    for gene in genes:
        localD = d[gene]
        for quadNuc in localD.keys():
            listOfDicts.append({'Hugo_Symbol': gene, 'quadNuc': quadNuc,
                                'totalNmut': localD[quadNuc]['totalNmut'], 'nOncogenicMut': localD[quadNuc]['nOncogenicMut']})
    return pd.DataFrame(listOfDicts)

#SUMMARIZES INFORMATION ABOUT THE SUSCEPTIBILITY OF DIFFERENT MOTIFS TO ONCOGENIC MUTATIONS
def summarize_data_about_quad_nuc_oncogenic_mut_susceptibility():
    validQuadNucs = set([start + change + end for start in ['A', 'C', 'T', 'G'] 
                        for change in ['CA', 'CG', 'CT', 'TA', 'TC', 'TG']
                        for end in ['A', 'C', 'T', 'G']])
    
    impact341OncogenicInfo = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/impact341SimulatedMutationDataSummary.tsv')
    impact341OncogenicInfo = impact341OncogenicInfo[impact341OncogenicInfo['quadNuc'].isin(validQuadNucs)]
    qNucDict = mutation_modeling_util.calculate_quadnuc_based_oncogenic_susceptibility_dict(impact341OncogenicInfo) 
    listOfDicts = []
    for key, value in qNucDict.items():
        listOfDicts.append({'QuadNuc': key, 'frac_oncogenic':value})
    df = pd.DataFrame(listOfDicts)
    df['orderingVal'] = df['QuadNuc'].apply(lambda x: 0 if x[1:3] == 'CA' else 1 if x[1:3] == 'CG' else 2 if x[1:3] == 'CT' else 3 if x[1:3] == 'TA' else 4 if x[1:3] == 'TC' else 5 if x[1:3] == 'TG' else None)
    df.to_csv('~/Desktop/WORK/dataForLocalPlotting/quadNucOncogenicInfo.tsv', sep='\t', index=False)
    

#summarize_data_about_observed_quad_nuc_oncogenic_mut_prevalence(mutationMaf)
#TELLS you the observed mutation prevalence of quad nucs
def summarize_data_about_observed_quad_nuc_oncogenic_mut_prevalence(mutDf):
    validQuadNucs = set([start + change + end for start in ['A', 'C', 'T', 'G'] 
                        for change in ['CA', 'CG', 'CT', 'TA', 'TC', 'TG']
                        for end in ['A', 'C', 'T', 'G']])
    oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
    listOfDicts = []
    impact341Genes = set(['ABL1', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BAP1', 'BARD1', 'BBC3', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CARD11', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79B', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSF1R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EP300', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERG', 'ESR1', 'ETV1', 'ETV6', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXP1', 'FUBP1', 'GATA1', 'GATA2', 'GATA3', 'GNA11', 'GNAQ', 'GNAS', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3B', 'HNF1A', 'HRAS', 'ICOSLG', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INPP4A', 'INPP4B', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAPK1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH6', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOR1', 'NF1', 'NF2', 'NFE2L2', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTRK1', 'NTRK2', 'NTRK3', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLE', 'PPP2R1A', 'PRDM1', 'PRKAR1A', 'PTCH1', 'PTEN', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAC1', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RUNX1', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SETD2', 'SF3B1', 'SH2D1A', 'SHQ1', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SOCS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SRC', 'STAG2', 'STK11', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TBX3', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP63', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'VHL', 'VTCN1', 'WT1', 'XIAP', 'XPO1', 'YAP1', 'YES1'])
    mutDf = mutDf[mutDf['Hugo_Symbol'].isin(impact341Genes)]
    for quadNuc in validQuadNucs:
        mutDataAtQuadNuc = mutDf[mutDf['quadNuc'] == quadNuc]
        nMutAtQuadNuc = mutDataAtQuadNuc.shape[0]
        nOncogneicMutAtQuadNuc = mutDataAtQuadNuc[mutDataAtQuadNuc['oncogenic'].isin(oncogenicMutColNames)].shape[0]
        fracObservedMutsOncogenic = 1.0*nOncogneicMutAtQuadNuc/nMutAtQuadNuc
        listOfDicts.append({'QuadNuc': quadNuc, 'frac_oncogenic': fracObservedMutsOncogenic, 'Nmut': nOncogneicMutAtQuadNuc})
    df = pd.DataFrame(listOfDicts)
    df['orderingVal'] = df['QuadNuc'].apply(lambda x: 0 if x[1:3] == 'CA' else 1 if x[1:3] == 'CG' else 2 if x[1:3] == 'CT' else 3 if x[1:3] == 'TA' else 4 if x[1:3] == 'TC' else 5 if x[1:3] == 'TG' else None)
    df.to_csv('~/Desktop/WORK/dataForLocalPlotting/quadNucOncogenicInfo_Observed.tsv', sep='\t', index=False)
    

##############TEMP NOAH
    
mutationMaf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/hypermutationAnalysisProjectData/hypermutationAnalysisProjectMaf.tsv')
mutationMaf['quadNuc'] = mutationMaf.apply(lambda row: mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)

possibleMutData = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/panImpactSimulatedMutationDataSummary.tsv')
oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
impact468 = set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'AMER1', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'GTF2I', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MLL4', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SETD8', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2']) 
for gene in impact468:
    geneInfo = possibleMutData[possibleMutData['Hugo_Symbol'] == gene]
    if geneInfo.shape[0] > 0:
        nOnc = sum(geneInfo['nOncogenicMut'])
        if nOnc == 0:
            nOncObs = mutationMaf[(mutationMaf['Hugo_Symbol'] == gene) & (mutationMaf['oncogenic'].isin(oncogenicMutColNames))].shape[0]
            print gene, ' NoncoObserved:', nOncObs,  ' CHROMOSOME:', geneInfo['Chromosome'].iloc[0]


print possibleMutData[possibleMutData['Hugo_Symbol'] == 'CIC']['Chromosome']


#ONLY DO ANALYSIS ON IMPACT 341 genes because these are the easiest to analyze 
impact341Genes = set(['ABL1', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BAP1', 'BARD1', 'BBC3', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CARD11', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79B', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSF1R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EP300', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERG', 'ESR1', 'ETV1', 'ETV6', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXP1', 'FUBP1', 'GATA1', 'GATA2', 'GATA3', 'GNA11', 'GNAQ', 'GNAS', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3B', 'HNF1A', 'HRAS', 'ICOSLG', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INPP4A', 'INPP4B', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAPK1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH6', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOR1', 'NF1', 'NF2', 'NFE2L2', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTRK1', 'NTRK2', 'NTRK3', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLE', 'PPP2R1A', 'PRDM1', 'PRKAR1A', 'PTCH1', 'PTEN', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAC1', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RUNX1', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SETD2', 'SF3B1', 'SH2D1A', 'SHQ1', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SOCS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SRC', 'STAG2', 'STK11', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TBX3', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP63', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'VHL', 'VTCN1', 'WT1', 'XIAP', 'XPO1', 'YAP1', 'YES1'])
mutationMaf341 = mutationMaf[mutationMaf['Hugo_Symbol'].isin(impact341Genes)]

impactSigs = pd.read_table(pathPrefix + '/ifs/res/taylorlab/impact_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
impactSigs['pid'] = impactSigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
impactSigs['cancer_type'] = impactSigs['pid'].apply(lambda x: cDict[x] if x in cDict else None)
poleEndometrials = set(impactSigs[(impactSigs['cancer_type'] == 'Endometrial Cancer') & (impactSigs['mean_10'] > .2) & (impactSigs['Nmut_Mb'] > 20)]['Tumor_Sample_Barcode'])
tmzGliomas = set(impactSigs[(impactSigs['cancer_type'] == 'Glioma') & (impactSigs['mean_11'] > .2) & (impactSigs['Nmut_Mb'] > 20)]['Tumor_Sample_Barcode'])
sig17High = set(impactSigs[(impactSigs['mean_17'] > .2) & (impactSigs['Nmut_Mb'] > 8)]['Tumor_Sample_Barcode'])

#DO COHORT SPECIFIC ANALYSES
cohortMaf1 = mutationMaf341[mutationMaf341.isin(poleEndometrials)] #todo make the cohort maf a value here
cohortMaf2 = mutationMaf341[mutationMaf341.isin(tmzGliomas)]
cohortMaf3 = mutationMaf341[mutationMaf341.isin(sig17High)]

d = mutation_modeling_util.initiate_gene_mut_mapping(simulatedMutData341)
df = convert_dict_to_df(d)

summarize_data_about_observed_quad_nuc_oncogenic_mut_prevalence(mutationMaf)
#FOR LOOP TO COMPARE OBSERVED VS EXPECTED

oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
oncogenicSusceptibilityDict = mutation_modeling_util.calculate_quadnuc_based_oncogenic_susceptibility_dict(possibleMutData)

cohortMaf = cohortMaf3
for case in set(cohortMaf['Tumor_Sample_Barcode']):
    caseSigs = impactSigs[impactSigs['Tumor_Sample_Barcode'] == case]
    #caseMaf = cohortMaf[co]
    if caseSigs.shape[0] > 0:
        signatureColumnNames = ['mean_' + str(i) for i in range(1,31)]
        signatureColsOnly = caseSigs[signatureColumnNames]
        renameDict = dict([('mean_' + str(i), 'Signature.' + str(i)) for i in range(1,31)])
        signatureColsOnly = signatureColsOnly.rename(columns=renameDict)
        
        decompositionDict = signatureColsOnly.to_dict(orient='records')[0]
        quadNucFractions = mutation_modeling_util.get_quadnuc_fracs_given_decomposition(decompositionDict, spectraPath = pathPrefix + '/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt')
        
        print mutation_modeling_util.get_expected_oncogenic_val_given_quadnuc_fractions(quadNucFractions, oncogenicSusceptibilityDict)


possibleMutData = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/impact341SimulatedMutationDataSummary.tsv')

mutation_modeling_util.get_chance_of_oncogenic_mutation_in_gene_at_quadnuc(possibleMutData, 'TCTG', 'APC')

for gene in set(possibleMutData['Hugo_Symbol']):
    print gene, sum(possibleMutData[possibleMutData['Hugo_Symbol'] == gene]['totalNmut'])

possibleMutData[possibleMutData['Hugo_Symbol'] == 'TERT']['quadNuc']

possibleMutData.columns.values

print sum(df['totalNmut']), sum(df['nOncogenicMut'])










































    

mafHere = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/annotatedOncoPlusHotspotMafAllImpact_trinuc')
#mafHere = pd.read_table('~/Desktop/OfflineStuffFor11_20Trip/annotatedOncoPlusHotspotMafAllImpact_trinuc')
#geneDist = pd.read_table('~/Desktop/OfflineStuffFor11_20Trip/impact_gene_trinuc_distributions.tsv')

mafHere['quadNuc'] = mafHere.apply(lambda row: mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)

#the background distribtuion of all observed mutations in impact
mafBackground = pd.read_table(pathPrefix + '/ifs/res/taylorlab/ang46/ext/dmp/mskimpact/mutation_data.txt')

#The background distribution of trinucs 
geneDistributionsDf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/impact_gene_trinuc_distributions.tsv')

#smallMafWithOncogenicAnnotations = mafHere.sample(3000, axis=0)
#smallMafBackgroundMaf = mafBackground.sample(3000, axis=0)
#smallMafWithOncogenicAnnotations['quadNuc'] = smallMafWithOncogenicAnnotations['quadNuc'] = smallMafWithOncogenicAnnotations.apply(lambda row: 
#   			mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)

    
#BIG LINE OF CODE THAT GETS ALL THE PROBABILITY MODELS WE NEED
#DO it quick with minimafs
#quadNucDict, quadNucDictAdjusted, geneOncogenicOrHotspotProbDict, geneMutProbInfo = mutation_modeling_util.initiate_models(smallMafWithOncogenicAnnotations, smallMafBackgroundMaf, geneDistributionsDf)
#Do it the long way with big ass full mafs
quadNucDict, quadNucDictAdjusted, geneOncogenicOrHotspotProbDict, geneMutProbInfo = mutation_modeling_util.initiate_models(mafHere, mafBackground, geneDistributionsDf)
mutation_modeling_util.write_model_information(quadNucDict=quadNucDict, quadNucDictAdjusted=quadNucDictAdjusted, geneOncogenicOrHotspotProbDict=geneOncogenicOrHotspotProbDict, geneMutProbInfo=geneMutProbInfo, write_dir='~/Desktop/mnt/ifs/work/taylorlab/friedman/myAdjustedDataFiles/mutationSimultationData',)

spectraD = mutationSigUtils.convert_spectrum_file_to_dict_of_dicts(spectrumFile=pathPrefix + '/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt')

#DO SIMULATION
vals = model_oncogenic_mut_acquisition(spectraD, quadNucDict, quadNucDictAdjusted, geneOncogenicOrHotspotProbDict, geneMutProbInfo, simulatedNmuts = 1000, simulatedKCases = 50, mutationMode='hotspot')

listOfDicts = []
#MAKE a summary dataframe using the long data format
for val in vals:
    localDict = dict()
    localDict['SigName'] = val[0]
    localDict['average'] = val[1]
    localDict['modelName'] = 'model1'
    listOfDicts.append(localDict)
    
    localDict = dict()
    localDict['SigName'] = val[0]
    localDict['average'] = val[2]
    localDict['modelName'] = 'model2'
    listOfDicts.append(localDict)
    
    localDict = dict()
    localDict['SigName'] = val[0]
    localDict['average'] = val[3]
    localDict['modelName'] = 'model3'
    listOfDicts.append(localDict)

dfData = pd.DataFrame(listOfDicts)

dfData.to_csv('~/Desktop/dataForLocalPlotting/simulationData.tsv', sep='\t', index=False)






quadNucDict, geneOncogenicOrHotspotProbDict = mutation_modeling_util.do_initiation(mafWithInfo, geneDistributions=pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/impact_gene_trinuc_distributions.tsv')


l = []
for i in range(100):
    l.append(mutation_modeling_util.pick_mutation_given_gene_and_quadNuc('ARID1A', 'TTCG', geneOncogenicOrHotspotProbDict, mode='hotspot'))
print Counter(l)


tmzSigs = gliomaSigs[gliomaSigs['mean_11'] > .25]
renameDict = {key:value for (key,value) in [('mean_' + str(i), 'Signature.' + str(i)) for i in range(1,31)]}
tmzSigs = tmzSigs.rename(columns=renameDict)
signatureCols = ['Signature.' + str(i) for i in range(1,31)]


val = asses_simulated_mutation_rates_by_signatures(sigMagnitudesDict, 100, spectraD, 1000, quadNucDict, geneOncogenicOrHotspotProbDict, mMode='oncogenic')
df = pd.DataFrame([])
df['signature'] = val.keys()
df['average'] = val.values()

df.to_csv('~/Desktop/noahQuickTest.tsv', sep='\t', index=False)

cases = set(tmzSigs['Tumor_Sample_Barcode'])
for case in cases:
    print mutDataDf[mutDataDf['Tumor_Sample_Barcode'] == case]['nMuts'].iloc[0], tmzSigs[tmzSigs['Tumor_Sample_Barcode'] == case]['Nmut'].iloc[0]
    #localD = tmzSigs[tmzSigs['Tumor_Sample_Barcode'] == case]
    #sigMagnitudesDict = localD[signatureCols].iloc[0].to_dict()
    
    #sigMagnitudesDict = {key:value for (key, value) in [('Signature.' + str(i), 0) for i in range(1,31)]}
    #sigMagnitudesDict['Signature.1'] = 1
    
    #mutation_modeling_util.do_k_simulations_of_mutations(sigMagnitudesDict, 10, spectraD, 1000, quadNucDict, geneOncogenicOrHotspotProbDict, mutMode='oncogenic')
    #spectraChosen, quadNucsChosen, genesChosen, mutationsChosen = mutation_modeling_util.pick_n_mutations_given_spectra(sigMagnitudesDict, spectraD, 1000, quadNucDict, geneOncogenicOrHotspotProbDict, mutMode='oncogenic')
    #print Counter(spectraChosen), Counter(quadNucsChosen), Counter(genesChosen), Counter(mutationsChosen)
    #break
    
    
#################################################################
#########Look at mutation landscapes across different mutation signatures
    
dfTest = pd.read_table(pathPrefix + '/ifs/res/taylorlab/ang46/ext/dmp/mskimpact/mutation_data.txt')


