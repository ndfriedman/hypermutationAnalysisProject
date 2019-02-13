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

cdsSizeInfo = analysis_utils.get_cds_size_targeted_by_impact(infoFilePath = pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/impact_gene_reference_signatures.tsv')
mutationMaf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/hypermutationAnalysisProjectData/hypermutationAnalysisProjectMaf.tsv')
cdsSizeInfo['PTPRT']

geneReferenceInfoDf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/impact_gene_reference_signatures.tsv')
geneReferenceInfoDf.columns.values['Hugo_Symbol']
oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])

mutationMaf['HGVSp_Short']
len(set(mutationMaf[(mutationMaf['Hugo_Symbol'] == 'PTEN') & (mutationMaf['oncogenic'].isin(oncogenicMutColNames)) & (mutationMaf['Variant_Type'] == 'SNP')]['Start_Position']))

mutationMaf['Variant_Type']


len(set(mutationMaf[(mutationMaf['Hugo_Symbol'] == 'JAK1') & (mutationMaf['oncogenic'].isin(oncogenicMutColNames)) & (mutationMaf['Variant_Type'] == 'SNP')]['Start_Position']))


listOfDicts = []
impactGenes = ['INPP4A', 'GNA11', 'MEF2B', 'FGF4', 'KEAP1', 'ESR1', 'MAPK1', 'NFKBIA', 'STAG2', 'NTHL1', 'TSC2', 'MGA', 'AGO2', 'CD79A', 'PIK3R2', 'KMT2B', 'ERF', 'HGF', 'PHOX2B', 'CCND1', 'GLI1', 'CDKN1B', 'RAB35', 'MLH1', 'BCL6', 'MSH2', 'MSH6', 'TNFAIP3', 'DUSP4', 'CXCR4', 'FLT3', 'INHBA', 'INHA', 'HIST1H3B', 'CDKN1A', 'SOX9', 'RRAS', 'TRAF2', 'RAC2', 'SMO', 'KNSTRN', 'MYOD1', 'FOXA1', 'RAF1', 'SESN2', 'LATS1', 'RARA', 'DNAJB1', 'H3F3B', 'KRAS', 'RRAS2', 'VHL', 'NOTCH2', 'PDGFRA', 'APC', 'MSI1', 'HNF1A', 'TBX3', 'CDK4', 'TMEM127', 'BRIP1', 'FGFR3', 'BARD1', 'CCND2', 'PALB2', 'CDH1', 'PDGFRB', 'FLT4', 'PIK3C3', 'SMAD2', 'RHEB', 'KMT2C', 'AXIN1', 'CREBBP', 'CCNE1', 'CDKN2C', 'PIK3R3', 'UPF1', 'MAP2K2', 'INPP4B', 'MAPK3', 'SMARCB1', 'EP300', 'EED', 'BRD4', 'NOTCH3', 'BIRC3', 'ACVR1', 'EPAS1', 'EPCAM', 'AKT3', 'KDR', 'PIK3CA', 'CTCF', 'CBL', 'CUL3', 'STAT3', 'DNMT3A', 'TP63', 'SDHA', 'MAP3K13', 'MSH3', 'RAD50', 'NBN', 'CDK6', 'PMS2', 'MAPKAP1', 'PIK3C2G', 'ERBB3', 'RB1', 'RAD51', 'IGF1R', 'ZFHX3', 'NCOR1', 'TP53', 'ERBB2', 'NUF2', 'EPHA5', 'PLK2', 'PIK3R1', 'RASA1', 'EGFR', 'NDUFB11', 'PRDM14', 'CDKN2B', 'NTRK2', 'NOTCH1', 'ATM', 'ARID5B', 'EIF4E', 'MYCN', 'FBXW7', 'CYSLTR2', 'FLT1', 'YAP1', 'MSI2', 'TCEB1', 'FLCN', 'ERCC3', 'CSF1R', 'GNAQ', 'PPARG', 'KIT', 'ERG', 'PREX2', 'BRAF', 'FANCC', 'PIK3CB', 'HIST1H2BD', 'HOXB13', 'U2AF1', 'FGFR4', 'DCUN1D1', 'STAT5B', 'SLX4', 'FGF19', 'REL', 'PRKCI', 'MST1R', 'NPM1', 'SOX17', 'RAD21', 'TSHR', 'INPPL1', 'TSC1', 'SPRED1', 'GREM1', 'RUNX1', 'ANKRD11', 'KMT2D', 'AXL', 'SDHAF2', 'CTLA4', 'INSR', 'IL7R', 'CDKN2A', 'IRS1', 'HIST1H3G', 'PPM1D', 'RPTOR', 'IGF1', 'AXIN2', 'MAP2K1', 'BCL2L1', 'ZRSR2', 'NUP93', 'BTK', 'EGFL7', 'TERT', 'HRAS', 'ERCC4', 'RPS6KB2', 'AURKA', 'YES1', 'CALR', 'GSK3B', 'PMAIP1', 'WHSC1L1', 'CD276', 'ABL1', 'FOXP1', 'ALOX12B', 'EZH2', 'POLE', 'FAM175A', 'PPP2R1A', 'SUZ12', 'MRE11A', 'EIF4A2', 'ARID1A', 'GTF2I', 'SOX2', 'PGR', 'SHQ1', 'TRAF7', 'STK11', 'CARM1', 'SMAD3', 'DNMT3B', 'CHEK2', 'HIST1H3I', 'IDH2', 'AMER1', 'FOXL2', 'SETD8', 'GRIN2A', 'IKZF1', 'HIST2H3D', 'PTCH1', 'PRKD1', 'SOCS1', 'WT1', 'BCL2', 'FGF3', 'RPS6KA4', 'ARID2', 'PDCD1', 'SF3B1', 'CENPA', 'LMO1', 'RAD51D', 'PNRC1', 'ASXL2', 'EPHA3', 'RAD51C', 'FIP1L1', 'MEN1', 'NF2', 'H3F3C', 'DNMT1', 'GATA2', 'SH2B3', 'PDPK1', 'JAK1', 'ERBB4', 'SMAD4', 'DICER1', 'HIST1H1C', 'CDC42', 'DROSHA', 'SMARCA4', 'TCF3', 'IDH1', 'STAT5A', 'ARID1B', 'GATA3', 'E2F3', 'SPOP', 'MALT1', 'AKT1', 'CTNNB1', 'ATR', 'PTPN11', 'MITF', 'PAK7', 'MAP2K4', 'TAP1', 'CRKL', 'NKX2-1', 'BLM', 'XIAP', 'RET', 'TNFRSF14', 'ERCC5', 'RAC1', 'PAK1', 'PTPRD', 'HIST1H3D', 'PPP4R2', 'PTPRS', 'RICTOR', 'HIST1H3A', 'BRCA1', 'NEGR1', 'PAX5', 'SRC', 'NF1', 'CASP8', 'FGFR2', 'RAD52', 'PRKAR1A', 'MAX', 'TGFBR2', 'PIK3CG', 'HIST1H3J', 'XRCC2', 'PLCG2', 'BABAM1', 'ELF3', 'SRSF2', 'HIST1H3E', 'WWTR1', 'NTRK3', 'TOP1', 'MTOR', 'CSF3R', 'SDCCAG8', 'FH', 'HIST3H3', 'PARP1', 'H3F3A', 'PARK2', 'IKBKE', 'MDM4', 'CDC73', 'RFWD2', 'IFNGR1', 'DDR2', 'SDHC', 'INSRR', 'RIT1', 'ROS1', 'FYN', 'MCL1', 'PRDM1', 'HIST1H3H', 'EPHA7', 'FAM46C', 'SHOC2', 'VTCN1', 'NRAS', 'SUFU', 'BCL10', 'PTP4A1', 'FUBP1', 'GNAS', 'SH2D1A', 'JUN', 'PTEN', 'RAD54L', 'NCOA3', 'BMPR1A', 'MUTYH', 'MPL', 'CCND3', 'RRAGC', 'STK40', 'PTPRT', 'ATRX', 'PIM1', 'PPP6C', 'TET1', 'MED12', 'DAXX', 'ID3', 'KLF4', 'AR', 'TAP2', 'TGFBR1', 'NOTCH4', 'STK19', 'KDM5C', 'SDHB', 'SDHD', 'ASXL1', 'SYK', 'SPEN', 'IRS2', 'MDC1', 'GATA1', 'HLA-A', 'ARAF', 'PIK3CD', 'ERRFI1', 'RBM10', 'DIS3', 'KDM6A', 'MYC', 'BCOR', 'FOXO1', 'EIF1AX', 'TET2', 'TEK', 'BRCA2', 'GPS2', 'NKX3-1', 'IRF4', 'CDK8', 'CRLF2', 'CD274', 'JAK2', 'TP53BP1', 'LATS2', 'WHSC1', 'SMYD3', 'ALK', 'FANCA', 'ERCC2', 'AKT2', 'CD79B', 'BCL2L11', 'PBRM1', 'SMARCD1', 'MYD88', 'ETV6', 'CARD11', 'NFE2L2', 'MYCL', 'PDCD1LG2', 'MET', 'EPHB1', 'CYLD', 'TMPRSS2', 'DOT1L', 'AC129492.6', 'MAP3K1', 'KDM5A', 'XPO1', 'SOS1', 'OBSL1', 'ETV1', 'FAM58A', 'ICOSLG', 'RNF43', 'SETD2', 'HLA-B', 'CBFB', 'RHOA', 'RECQL', 'IL10', 'FGFR1', 'RECQL4', 'EZH1', 'CHEK1', 'IGF2', 'SESN1', 'CSDE1', 'NSD1', 'POLD1', 'PMS1', 'FAT1', 'HIST1H3F', 'CDK12', 'BBC3', 'RP11-211G3.3', 'MST1', 'JAK3', 'BAP1', 'MDM2', 'RYBP', 'RXRA', 'RAD51B', 'CEBPA', 'RTEL1', 'LYN', 'VEGFA', 'NTRK1', 'KMT2A', 'SESN3', 'HIST1H3C', 'TIMM8B', 'TCF7L2', 'B2M', 'CIC', 'AURKB', 'MFSD11']
for gene in impactGenes:
    nDistinctOncogenicSNPs = len(set(mutationMaf[(mutationMaf['Hugo_Symbol'] == gene) & (mutationMaf['oncogenic'].isin(oncogenicMutColNames)) & (mutationMaf['Variant_Type'] == 'SNP')]['Start_Position']))
    cdsLength = cdsSizeInfo[gene]
    listOfDicts.append({'gene': gene, 'cdsLength': cdsLength, 'nDistinctOncogenicSnps': nDistinctOncogenicSNPs})
 
mutationMaf[mutationMaf['Hugo_Symbol'] == 'KMT2C']['Start_Position']    
    
len(set(mutationMaf[(mutationMaf['Hugo_Symbol'] == 'PTEN')&(mutationMaf['oncogenic'].isin(oncogenicMutColNames))&(mutationMaf['Variant_Classification'] == 'Missense_Mutation')]['HGVSp_Short']))


mutationMaf.columns.values

df = pd.DataFrame(listOfDicts)
df['orderingVal'] = df['nDistinctOncogenicSnps']
df['logCdsValue'] = df['cdsLength'].apply(lambda x: math.log(x, 2))

df = df[df['orderingVal'].notnull()]

df.to_csv('~/Desktop/WORK/dataForLocalPlotting/geneOncogenicInfo.tsv', sep='\t', index=False)




mutationSigUtils.create_reference_four_nuc('GTA', 'C', 'T', 'SNP')




