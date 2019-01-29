#written by Noah Friedman (a template for scripts to be excuted in the spyder environment
#this script is used to create simple fisher test summary statistics for volcano plots of gene mutation prevalence across different poulations
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
import scipy
import math
import double_mutation_analysis_util

def fisher_test_for_gene(df, gene, sepCol='isHypermutator', mutCategoryColumn = 'mutKnowOncogenicOrTruncating'):
    
    df = df[df[mutCategoryColumn] == True]
    
    yesDf = df[df[sepCol] == True]
    noDf = df[df[sepCol] != True]
    nYes = len(set(yesDf['Tumor_Sample_Barcode']))
    nNo = len(set(noDf['Tumor_Sample_Barcode']))
    
    nGeneYes = len(set(yesDf[yesDf['Hugo_Symbol'] == gene]['Tumor_Sample_Barcode']))
    nNotGeneYes = nYes - nGeneYes
    
    nGeneNo = len(set(noDf[noDf['Hugo_Symbol'] == gene]['Tumor_Sample_Barcode']))
    nNotGeneNo = nNo - nGeneNo
    
    #Column 1 of table is hypermutator, Column 2 is not hypermutators.  Row 1 is gene is mutated, Row 2 gene is not mutated
    table = [[nGeneYes, nGeneNo], [nNotGeneYes, nNotGeneNo]]
    
    #ITS POSSIBLE THE TABLE IS WRONG!!!
    oddsratio, pvalue = scipy.stats.fisher_exact(table)
    if oddsratio == 0 or pvalue == 0:
        return None, None
    return math.log(oddsratio, 2), -1*math.log(pvalue, 10)

#DOES THIS STUFF IN MODE 2 which is working on a doublet formatted Df
def fisher_test_for_gene_mode_2(df, geneCol, sepCol='isHypermutator'):

    yesDf = df[df[sepCol] == True]
    noDf = df[df[sepCol] != True]
    nYes = len(set(yesDf['Tumor_Sample_Barcode']))
    nNo = len(set(noDf['Tumor_Sample_Barcode']))
    
    nGeneYes = len(set(yesDf[yesDf[geneCol] == True]['Tumor_Sample_Barcode']))
    nNotGeneYes = nYes - nGeneYes
    
    nGeneNo = len(set(noDf[noDf[geneCol] == False]['Tumor_Sample_Barcode']))
    nNotGeneNo = nNo - nGeneNo
    
    #Column 1 of table is hypermutator, Column 2 is not hypermutators.  Row 1 is gene is mutated, Row 2 gene is not mutated
    table = [[nGeneYes, nGeneNo], [nNotGeneYes, nNotGeneNo]]
    
    print geneCol, table
    
    #ITS POSSIBLE THE TABLE IS WRONG!!!
    oddsratio, pvalue = scipy.stats.fisher_exact(table)
    if oddsratio == 0 or pvalue == 0:
        return None, None
    return math.log(oddsratio, 2), -1*math.log(pvalue, 10)




#############################################
hypermuationThreshold = 80

maf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/hypermutationAnalysisProjectData/hypermutationAnalysisProjectMaf.tsv')
sigs = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/hypermutationAnalysisProjectData/hypermutationAnalysisProjectSignatures.tsv')
#sigs = mutationSigUtils.merge_signature_columns(sigs)
#sigs['dominantSignature'] = sigs.apply(lambda row: mutationSigUtils.get_dominant_signature(row.to_dict(), cols=None), axis=1)
  
sigs['pid'] = sigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
sigs['cancer_type'] = sigs['pid'].apply(lambda x: cDict[x] if x in cDict else None)

hypermutatorIds = set(sigs[sigs['Nmut_Mb'] > hypermuationThreshold]['Tumor_Sample_Barcode'])
endometrialIds = set(sigs[sigs['cancer_type'] == 'Endometrial Cancer']['Tumor_Sample_Barcode'])
endometrialHypermutators = hypermutatorIds & endometrialIds

phasedDf = pd.read_table(pathPrefix + '/home/ang46/luna/projects/multiplets/data/processed_data/data_phased_mutationphaser_v0.3.1_annotated.txt')
phasedDf['isHypermutator'] = phasedDf['Tumor_Sample_Barcode'].apply(lambda x: True if x in hypermutatorIds else False)
#phasedDfEndometrialHypers = phasedDf[phasedDf['Tumor_Sample_Barcode'].isin(endometrialHypermutators)]
endometrialPhasedDf = phasedDf[phasedDf['Tumor_Sample_Barcode'].isin(endometrialIds)]

endometrialPhasedDf[(endometrialPhasedDf['Hugo_Symbol'] == 'PTEN') & (endometrialPhasedDf['isHypermutator'] == False)]

endometrialMaf = maf[maf['Tumor_Sample_Barcode'].isin(endometrialIds)]

endometrialNormalMaf = maf[maf['Tumor_Sample_Barcode'].isin(endometrialIds - hypermutatorIds)]
endometrialHyperMaf = maf[maf['Tumor_Sample_Barcode'].isin(endometrialIds & hypermutatorIds)]

len(set(endometrialHyperMaf['Tumor_Sample_Barcode']))

doubCntr = 0
for tid in set(endometrialHyperMaf['Tumor_Sample_Barcode']):
    caseMaf = maf[maf['Tumor_Sample_Barcode'] == tid]
    if caseMaf[caseMaf['Hugo_Symbol'] == 'PTEN'].shape[0] > 1:
        doubCntr +=1
print doubCntr

table = [[32, 151], [9, 558]]
table= [[29, 198], [12, 141]]
    #ITS POSSIBLE THE TABLE IS WRONG!!!
oddsratio, pvalue = scipy.stats.fisher_exact(table)
math.log(oddsratio, 2), -1*math.log(pvalue, 10)

returnedDf = double_mutation_analysis_util.create_double_mutation_summary_vanilla_maf(endometrialHyperMaf)

Counter(returnedDf['CTNNB1_doubleOncogenic'])
Counter(returnedDf['CTNNB1_double'])

dfGenes = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/impact_gene_reference_signatures.tsv')
dfGenes.columns.values   
set(dfGenes['Hugo_Symbol'])
print sum(list(dfGenes[(dfGenes['Hugo_Symbol'] == 'PTEN')].iloc[0])[3:])
 
print Counter(dict(zip(dfGenes['Hugo_Symbol'], dfGenes['cds_length']))).most_common(600)

dfProbeDesign = pd.read_table(pathPrefix + '/ifs/depot/resources/targets/IMPACT410_hg19/IMPACT410_RC2.berger')
print sum(list(dfProbeDesign[(dfProbeDesign['GeneID'] == 'PTEN')  & (dfProbeDesign['Category'] == 'Exon')]['Length']))

dfProbeDesign.columns.values

#########################################3

len(set(endometrialMaf['Tumor_Sample_Barcode']))
print maf[maf['Hugo_Symbol'] == 'PTEN']

doubleMutSummaryDf = double_mutation_analysis_util.create_double_mutation_summary_alex_phasing(endometrialPhasedDf)

impactGenes = set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KMT2B', 'KMT5A', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51C', 'RAD51L1', 'RAD51L3', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2'])

listOfDicts = []
for gene in impactGenes:
    log2oddsratio, log10pvalue = fisher_test_for_gene_mode_2(doubleMutSummaryDf, geneCol=gene + '_double', sepCol='isHypermutator')
    #print gene, log2oddsratio, log10pvalue
    #log2oddsratio, log10pvalue = fisher_test_for_gene(maf, gene)
    listOfDicts.append({'Gene': gene, 'log2oddsratio': log2oddsratio, 'log10pvalue': log10pvalue})

dfDoubleWithOncogenic = pd.DataFrame(listOfDicts)


for tid in set(phasedDfEndometrialHypers['Tumor_Sample_Barcode']):
   caseDf = phasedDfEndometrialHypers[phasedDfEndometrialHypers['Tumor_Sample_Barcode'] == tid]
   initShape = caseDf[caseDf['Hugo_Symbol'] == 'PTEN'].shape[0]
   if initShape == 0:
       print tid


print maf.columns.values

print len(ids)

endometrialHyperSigs = sigs[(sigs['cancer_type'] == 'Endometrial Cancer') & (sigs['Nmut_Mb'] > hypermuationThreshold)]

poleHyperMuts['quadNuc'] = poleHyperMuts.apply(lambda row: 
   			mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)


poleHyperMuts.columns.values
double









#maf['pid'] = maf['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
#cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
#maf['cancer_type'] = maf['pid'].apply(lambda x: cDict[x] if x in cDict else None)

print sigs[sigs['Tumor_Sample_Barcode'] == '']

gliomaTMZIds = sigs[(sigs['cancer_type'] == 'Glioma') & (sigs['Nmut_Mb'] > 80)]['Tumor_Sample_Barcode']

for i in v:
    print i

hypermutationThreshold = 80
maf['isHypermutator'] = maf['Nmut_Mb'].apply(lambda x: True if x >= hypermutationThreshold else False)

maf['is-a-hotspot'] = maf['is-a-hotspot'].apply(lambda x: True if x == 'Y' else False)

impactGenes = set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KMT2B', 'KMT5A', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51C', 'RAD51L1', 'RAD51L3', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2'])

listOfDicts = []

for gene in impactGenes:
    print gene
    log2oddsratio, log10pvalue = fisher_test_for_gene(maf, gene, mutCategoryColumn='is-a-hotspot')
    #log2oddsratio, log10pvalue = fisher_test_for_gene(maf, gene)
    listOfDicts.append({'Gene': gene, 'log2oddsratio': log2oddsratio, 'log10pvalue': log10pvalue})

df = pd.DataFrame(listOfDicts)
geneLengthDict = analysis_utils.get_gene_length_info(bedFilePath = pathPrefix + '/ifs/res/pwg/data/gencode/gencode.v19.all_gene_bounds.bed')

print np.nanmedian(geneLengthDict.values())

geneLengthDict['BRCA2']

print Counter(geneLengthDict).most_common(500)
max(geneLengthDict.iteritems(), key=operator.itemgetter(1))

df['log2GeneLength'] = df['Gene'].apply(lambda x: math.log(geneLengthDict[x],2) if x in geneLengthDict else None)
#df['label'] = df['Gene']
df['label'] = None
genesToDisplay = set(['KRAS', 'EGFR', 'POLE', 'MLH1', 'MSH2', 'MSH6', 'PMS2', 'TP53', 'PIK3CA', 'NF1', 'ARID1A', 'BRCA1', 'BRCA2', 'MUTYH', 'APC', 'ATM'])
#genesToDisplay = set(['KRAS', 'EGFR', 'POLE', 'MLH1', 'MSH2', 'MSH6', 'PMS2', 'TP53', 'PIK3CA', 'NF1', 'ARID1A', 'BRCA1', 'BRCA2', 'MUTYH', 'APC', 'ATM', 'TP63', 'STK19'])
df['label'] = df['label'].apply(lambda x: x if x in genesToDisplay else None)

df['label'] = df.apply(lambda row: row['Gene'] if row['log2oddsratio'] > 3 and row['log10pvalue'] > 6 or row['log2oddsratio'] < -.5 else row['label'], axis=1)

df.to_csv('~/Desktop/dataForLocalPlotting/volcanoPlotData.tsv', sep='\t', index=False)

df[df['log2oddsratio'] < 0]

print df.shape

#df[df['log2oddsratio'] > 4.5]
df[df['log10pvalue'] > 50]




######################################################33

def make_summary_dict(df, title, cancerTypesToAnalyze = set(['Non-Small Cell Lung Cancer', 'Melanoma', 'Colorectal Cancer', 'Bladder Cancer', 'Breast Cancer', 'Endometrial Cancer', 'Skin Cancer, Non-Melanoma', 'Soft Tissue Sarcoma', 'Glioma'])):
    dfSize = df.shape[0]
    localDict = {}
    countSum = 0
    for cancerType in cancerTypesToAnalyze:
        curCount = df[df['cancer_type'] == cancerType].shape[0]
        localDict[cancerType] = curCount
        countSum += curCount
    localDict['other'] = dfSize - countSum
    localDict['range'] = title
    return localDict



























#TODO MOVE THIS INTO ANOTHER SCRIPT

cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/dmp/mskimpact/data_clinical_sample.txt')

sigs['pid'] = sigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
sigs['cancer_type'] = sigs['pid'].apply(lambda x: cDict[x] if x in cDict else None)

sigs.columns.values

hypermutatorSigs = sigs[sigs['Nmut_Mb'] > 25]
hypermutatorSigs = hypermutatorSigs[hypermutatorSigs['cancer_type'].notnull()]

cancerTypesToAnalyze = set(['Non-Small Cell Lung Cancer', 'Melanoma', 'Colorectal Cancer', 'Bladder Cancer', 'Breast Cancer', 'Endometrial Cancer', 'Skin Cancer, Non-Melanoma', 'Soft Tissue Sarcoma', 'Glioma'])

cohort1 = hypermutatorSigs[(hypermutatorSigs['Nmut_Mb'] > 25) & (hypermutatorSigs['Nmut_Mb'] <= 50)]
cohort2 = hypermutatorSigs[(hypermutatorSigs['Nmut_Mb'] > 50) & (hypermutatorSigs['Nmut_Mb'] <= 75)]
cohort3 = hypermutatorSigs[(hypermutatorSigs['Nmut_Mb'] > 75) & (hypermutatorSigs['Nmut_Mb'] <= 100)]
cohort4 = hypermutatorSigs[(hypermutatorSigs['Nmut_Mb'] > 100) & (hypermutatorSigs['Nmut_Mb'] <= 200)]
cohort5 = hypermutatorSigs[(hypermutatorSigs['Nmut_Mb'] > 200) & (hypermutatorSigs['Nmut_Mb'] <= 400)]
cohort6 = hypermutatorSigs[hypermutatorSigs['Nmut_Mb'] > 400]

listOfDictsHere = [make_summary_dict(cohort1, "25-50mutMb"),
                   make_summary_dict(cohort2, "50-75mutMb"),
                   make_summary_dict(cohort3, "75-100mutMb"),
                   make_summary_dict(cohort4, "100-200mutMb"),
                   make_summary_dict(cohort5, "200-400mutMb"),
                   make_summary_dict(cohort6, "400+mutMb")]

df = pd.DataFrame(listOfDictsHere)
dfMelted = pd.melt(df, id_vars=['range'])

dfMelted.to_csv('~/Desktop/dataForLocalPlotting/cohortCancerTypeSummary.tsv', sep='\t', index=False)



for v in list(sigs[(sigs['Nmut_Mb'] > 100) & (sigs['mean_10'] > .25)]['Tumor_Sample_Barcode']):
    print v




print cohort1.shape, cohort2.shape, cohort3.shape, cohort4.shape, cohort5.shape, cohort6.shape

















