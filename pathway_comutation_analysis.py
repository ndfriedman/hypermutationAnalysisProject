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

#NOTE ONLY USING IMPACT PANEL GENES

pathways = {
        'RTK-RAS': ['EGFR', 'ERBB2', 'ERBB3', 'ERBB4', 'MET', 'PDGFRA', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'KIT', 'IGFR1', 
                    'RET', 'ROS1', 'ALK', 'FLT3', 'JAK2', 'NTRK1-3', 'CBL', 'ERRfI1', 'ABL1', 'SOS1', 'NF1', 'RASA1', 'PTPN11',
                    'KRAS', 'HRAS', 'NRAS', 'RIT1', 'ARAF', 'BRAF', 'RAF1', 'RAC1', 'MAPK1', 'MAP2K1', 'MAP2K1', 'MAP2K2'],
                    
        'PI3K': ['PTEN', 'INPP4B', 'PIK3CA', 'PIK3CB', 'PIK3R2', 'PIK3R1', 'PIK3R3', 'AKT1', 'AKT2', 'AKT3', 'PPP2R1A',
                 'STK11', 'TSC1', 'TSC2', 'RHEB', 'RPTOR', 'MTOR', 'RICTOR'],
                 
        'P53': ['CDKN2A', 'MDM2', 'MDM4', 'TP53', 'ATM', 'CHEK2', 'RPS6KA3'],
        
        'CELL_CYCLE': ['CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CCNE1', 'RB1', 'E2F1', 'E2F3', 'CCND1', 'CCND2', 'CCND3', 'CDK2', 'CDK4', 'CDK6'],
        
        'MYC': ['MYC', 'MYCN', 'MYCL', 'MAX', 'MGA', 'MNT', 'MIXI1', 'MXD1', 'MXD3', 'MXD4', 'MLX', 'MLXIP', 'MLXIPL'],
        
        'WNT': ['WIF1', 'SFRP1', 'RNF43', 'ZNRF3', 'GSK3B', 'CTNNB1', 'APC', 'AXIN1', 'AXIN2', 'AMER1', 'TCF7', 'TCF7L2'],
        
        'HIPPO': ['FAT1', 'LATS1', 'LATS2', 'NF2', 'YAP1'],
        
        'NOTCH': ['NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'FBXW7', 'CREBBP', 'EP300', 'NCOR1', 'NCOR2', 'KDM5A'],
        
        'TGFB': ['TGFBR1', 'TGFBR2', 'SMAD2', 'SMAD3', 'SMAD4', 'ACVR1'],
        
        'NRF1': ['KEAP1' , 'CUL3', 'NFE2L2'],
        
        'PI3KA-ENDOMETRIAL': ['PIK3CA', 'PTEN', 'PIK3R1'],
        
        'RTK-ENDOMETRIAL': ['FGFR2', 'ERBB2', 'KRAS', 'SOX17', 'GSK3B', 'FBXW7', 'CTNNB1'],
        
        'P53-COLORECTAL': ['ATM', 'TP53'],
        
        'WNT-COLORECTAL': ['CTNNB1', 'APC', 'AXIN2', 'TCF7L2', 'SOX9', 'FBXW7', 'ARID1A', 'MYC'],
        
        'RTK/RAS-COLORECTAL': ['ERBB2', 'ERBB3', 'NRAS', 'KRAS', 'BRAF'],
        
        'PTEN/PIK3CA-COLORECTAL': ['IGF2', 'IGF1R', 'PTEN', 'IRS2', 'PIK3CA', 'PIK3R1'],
        
        'P53-GLIOMA': ['MDM2', 'MDM4', 'TP53'],
        
        'RTK/RAS-GLIOMA': ['EGFR', 'PDGFRA', 'MET', 'FGFR', 'PTEN', 'RAS', 'NF1', 'BRAF'],
        
        'CELL_CYCLE-GLIOMA': ['CDKN2A', 'CDKN2B', 'CDKN2C', 'CDK4', 'CDK6', 'RB1']
        
            }

def create_mut_not_mut_matrix(maf):
    totalCountMat = np.array([])
    genes = ['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KMT2B', 'KMT5A', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'KMT2A', 'KMT2B', 'KMT2C', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51C', 'RAD51L1', 'RAD51L3', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2']
    for case in set(maf['Tumor_Sample_Barcode']): 
        caseMaf = maf[maf['Tumor_Sample_Barcode'] == case]
        mutatedGenes = set(caseMaf['Hugo_Symbol'])
        vec = np.array([1 if x in mutatedGenes else 0 for x in genes])[np.newaxis] #create a numpy matrix of 0s and 1s for occurence/absence
        coOccurenceMatrix = vec*vec.T #do matrix multiplication to get a matrix of mutation co-occurence
        
        if totalCountMat.shape == (0,): #if its the very first version 
            totalCountMat = coOccurenceMatrix
        else:
            totalCountMat += coOccurenceMatrix
    #now return a data frame that we can write and use in r
    #add a name column based on the indicies
    df = pd.DataFrame(data=totalCountMat, index=genes, columns=genes)
    df['names'] = df.index
    df = df.drop_duplicates(subset = ['names'])
    return df

def get_matrix_for_pathway(df, pathwayName, genes, pathways):
    dfPathway = df[df['names'].isin(pathways[pathwayName])]
    dfPathway = dfPathway[[i for i in pathways[pathwayName] if i in dfPathway.columns.values] + ['names']]
    
    #DO dataframe manipulation so it will plot nicely in R
    dfPathway = dfPathway.melt(id_vars=['names'])
    dfPathway['variable'] = dfPathway.apply(lambda row: row['variable'] if genes.index(row['variable']) >= genes.index(row['names']) else None, axis=1)
    dfPathway = dfPathway[dfPathway['variable'].notnull()]
    return dfPathway



maf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/hypermutationAnalysisProjectData/hypermutationAnalysisProjectMaf.tsv')
sigs = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/hypermutationAnalysisProjectData/hypermutationAnalysisProjectSignatures.tsv')
sigs['pid'] = sigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
sigs['cancer_type'] = sigs['pid'].apply(lambda x: cDict[x] if x in cDict else None)

genes = ['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KMT2B', 'KMT5A', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'KMT2A', 'KMT2B', 'KMT2C', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51C', 'RAD51L1', 'RAD51L3', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2']


endometrialHyperIds = set(sigs[(sigs['cancer_type'] == 'Endometrial Cancer') & (sigs['Nmut_Mb'] > 50)]['Tumor_Sample_Barcode'])
endometrialNotHyperIds = set(sigs[(sigs['cancer_type'] == 'Endometrial Cancer') & (sigs['Nmut_Mb'] < 50)]['Tumor_Sample_Barcode'])
oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
oncogenicMuts = maf[maf['oncogenic'].isin(oncogenicMutColNames)]

#endometrialOncogenicMaf = oncogenicMuts[oncogenicMuts['Tumor_Sample_Barcode'].isin(endometrialHyperIds)]
endometrialOncogenicMaf = oncogenicMuts[oncogenicMuts['Tumor_Sample_Barcode'].isin(endometrialNotHyperIds)]

coOccurenceDfEndometrial = create_mut_not_mut_matrix(endometrialOncogenicMaf)
dfPI3KEndometrial = get_matrix_for_pathway(coOccurenceDfEndometrial, 'PI3K', genes, pathways)
dfPI3KEndometrial.to_csv('~/Desktop/WORK/dataForLocalPlotting/coOccurenceData.tsv', index = False, sep='\t')

colorectalHyperIds = set(sigs[(sigs['cancer_type'] == 'Colorectal Cancer') & (sigs['Nmut_Mb'] > 100)]['Tumor_Sample_Barcode'])
colorectalNotHyperIds = set(sigs[(sigs['cancer_type'] == 'Colorectal Cancer') & (sigs['Nmut_Mb'] < 50)]['Tumor_Sample_Barcode'])
colorectalOncogenicMaf = oncogenicMuts[oncogenicMuts['Tumor_Sample_Barcode'].isin(colorectalHyperIds)]
colorectalOncogenicMaf = oncogenicMuts[oncogenicMuts['Tumor_Sample_Barcode'].isin(colorectalNotHyperIds)]
colorectalOccurenceDf = create_mut_not_mut_matrix(colorectalOncogenicMaf)

dfWntColon = get_matrix_for_pathway(colorectalOccurenceDf, 'WNT-COLORECTAL', genes, pathways)
dfWntColon.to_csv('~/Desktop/WORK/dataForLocalPlotting/coOccurenceData.tsv', index = False, sep='\t')


gliomaHyperIds = set(sigs[(sigs['cancer_type'] == 'Glioma') & (sigs['Nmut_Mb'] > 50)]['Tumor_Sample_Barcode'])
gliomaOncogenicMaf = oncogenicMuts[oncogenicMuts['Tumor_Sample_Barcode'].isin(gliomaHyperIds)]
coOccurenceDfGlioma = create_mut_not_mut_matrix(gliomaOncogenicMaf)
dfRTKRASGlioma = get_matrix_for_pathway(coOccurenceDfGlioma, 'RTK-RAS', genes, pathways)
dfRTKRASGlioma.to_csv('~/Desktop/WORK/dataForLocalPlotting/coOccurenceData.tsv', index = False, sep='\t')




endometrialOncogenicMaf[endometrialOncogenicMaf['Hugo_Symbol'] == 'BRAF']
genes.index('APC')

dfPI3KEndometrial[dfPI3KEndometrial['names'] == 'BRAF']

