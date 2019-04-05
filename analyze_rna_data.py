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

dfRna = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/rnaData/rnaDataTableFromTCGA')
len(dfRna.columns.values)

#######LOAD IN SIGNATURE INFO FILE BY FILE

endometrialTCGASigs = pd.read_table(pathPrefix + '/ifs/res/taylorlab/jonssonp/msk_impact_gml_som/tcga/mutsig/UCEC.mc3.v0.2.8.PUBLIC.maf.2018-01-29.oncokb.vep.mutsig.txt')
endometrialTCGASigs['tcgaIdShort'] = endometrialTCGASigs['Sample Name'].apply(lambda x: x[:16])
endometrialIds = set(endometrialTCGASigs['tcgaIdShort'])

gliomaTCGASigs = pd.read_table(pathPrefix + '/ifs/res/taylorlab/jonssonp/msk_impact_gml_som/tcga/mutsig/GBM.mc3.v0.2.8.PUBLIC.maf.2018-01-29.oncokb.vep.mutsig.txt')
gliomaTCGASigs['tcgaIdShort'] = gliomaTCGASigs['Sample Name'].apply(lambda x: x[:16])
gliomaIds = set(gliomaTCGASigs['tcgaIdShort'])

breastTCGASigs = pd.read_table(pathPrefix + '/ifs/res/taylorlab/jonssonp/msk_impact_gml_som/tcga/mutsig/BRCA.mc3.v0.2.8.PUBLIC.maf.2018-01-29.oncokb.vep.mutsig.txt')
breastTCGASigs['tcgaIdShort'] = breastTCGASigs['Sample Name'].apply(lambda x: x[:16])
breastIds = set(breastTCGASigs['tcgaIdShort'])

pancreaticTCGASigs = pd.read_table(pathPrefix + '/ifs/res/taylorlab/jonssonp/msk_impact_gml_som/tcga/mutsig/PAAD.mc3.v0.2.8.PUBLIC.maf.2018-01-29.oncokb.vep.mutsig.txt')
pancreaticTCGASigs['tcgaIdShort'] = pancreaticTCGASigs['Sample Name'].apply(lambda x: x[:16])
pancreaticIds = set(pancreaticTCGASigs['tcgaIdShort'])

arbitraryHypermutationThreshold = 3000
endometrialHyperIds = set(endometrialTCGASigs[endometrialTCGASigs['Number of Mutations'] > arbitraryHypermutationThreshold]['Sample_ID'])



#Do a bit of babadookery to match ids and genes
idsRNA = set([x[:12] for x in dfRna.columns.values])
matchedIdStrings = idsRNA & endometrialHyperIds
rnaCols = [x for x in dfRna.columns.values if x[:12] in matchedIdStrings] 
dfRna['Hugo_Symbol'] = dfRna['gene_id'].apply(lambda x: x.split('|')[0])

#transpose it to work with it better
dfRna.index = dfRna['Hugo_Symbol']
dfRnaT = dfRna.drop(['Hugo_Symbol'], axis=1).transpose() #drop the hugo symbol column and transpose so our index(hugo symbol) is now the columns
dfRnaT['tcgaId'] = dfRnaT.index
dfRnaT = dfRnaT.drop(['gene_id']) #drop the first row which is useless stuff
dfRnaT['tcgaIdShort'] = dfRnaT['tcgaId'].apply(lambda x: x[:16])

impact468genes = set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'AMER1', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'GTF2I', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MLL4', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SETD8', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2']) 
genes = ['TP53', 'NF1', 'CTNNB1', 'PIK3CA', 'PTEN', 'KRAS', 'GAPDH', 'SMAD4', 'POLE']

#we need pd.to_numeric because the columns by default are of type object
#todo do this differently please
dfCases = dfRnaT[dfRnaT['tcgaId'].isin(set(rnaCols))]
for gene in impact468genes:
    if gene in dfCases.columns.values:
        print gene, np.nanmean(pd.to_numeric(dfCases[gene]))    
    
    np.nanmedian(pd.to_numeric(dfRnaT['PTEN']))

print np.nanmedian(pd.to_numeric(dfRnaT[dfRnaT['tcgaIdShort'].isin(breastIds)]['KRAS']))
print np.nanmedian(pd.to_numeric(dfRnaT[dfRnaT['tcgaIdShort'].isin(pancreaticIds)]['KRAS']))

endometrialIds

print dfRnaT[dfRnaT['PTEN'].notnull()].shape

print dfRnaT['BAF250'].dtype

"""for col in rnaCols:
    print '\n'
    print col
    geneRna = set(dfRna['Hugo_Symbol'])
    for gene in genes:
        if gene in geneRna:
            print gene, dfRna[dfRna['Hugo_Symbol'] == gene][col].iloc[0]"""
        
#PLOT THE OVERALL DISTRIBUTION OF RNA SEQ LEVELS
#where is NF1 along that distribution in endo vs not
#box plots comparing nf1 in endometrial bialleics and others
            
        
#print dfRna['TCGA-A5-A2K5-01A-11R-A180-07']
#print dfRna.columns.values
print np.nanmean(dfRna['TCGA-A5-A2K5-01A-11R-A180-07'])