#written by Noah Friedman (a template for scripts to be excuted in the spyder environment
import sys
import argparse
import os
import pandas as pd
import numpy as np
import numpy

from collections import Counter

pathPrefix = ''
if os.getcwd() == '/Users/friedman/Desktop/mnt':
	pathPrefix = '/Users/friedman/Desktop/mnt'

sys.path.append(pathPrefix + '/ifs/work/taylorlab/friedman/myUtils')
import analysis_utils 
import mutationSigUtils 
import maf_analysis_utils

mafAnnoMaf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/All.dmp_somatic_data_mutations_unfiltered.mafAnno.oncokb.hotspots.maf', sep='\t')
mafAnnoMaf['vaf'] = mafAnnoMaf.apply(lambda row: 1.0*row['t_alt_count']/(row['t_alt_count'] + row['t_ref_count']) if row['t_alt_count'] + row['t_ref_count'] > 0 else None, axis = 1)

impactSigs = pd.read_table(pathPrefix + '/ifs/res/taylorlab/impact_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
impactSigs['pid'] = impactSigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
impactSigs['cancer_type'] = impactSigs['pid'].apply(lambda x: cDict[x] if x in cDict else None)

hypermutationThresh = 50
endometrialHyperIds = set(impactSigs[(impactSigs['cancer_type'] == 'Endometrial Cancer') & (impactSigs['Nmut_Mb'] > hypermutationThresh)]['Tumor_Sample_Barcode'])
gliomaHyperIds = set(impactSigs[(impactSigs['cancer_type'] == 'Glioma') & (impactSigs['Nmut_Mb'] > hypermutationThresh)]['Tumor_Sample_Barcode'])
colorectalHyperIds = set(impactSigs[(impactSigs['cancer_type'] == 'Colorectal Cancer') & (impactSigs['Nmut_Mb'] > hypermutationThresh)]['Tumor_Sample_Barcode'])

genes = set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KMT2B', 'KMT5A', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51C', 'RAD51L1', 'RAD51L3', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2'])
oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])

#cohortMaf = mafAnnoMaf[mafAnnoMaf['Tumor_Sample_Barcode'].isin(endometrialHyperIds) & (mafAnnoMaf['oncogenic'].isin(oncogenicMutColNames))]
cohortMaf = mafAnnoMaf[mafAnnoMaf['Tumor_Sample_Barcode'].isin(endometrialHyperIds)]
cohortMaf = mafAnnoMaf[mafAnnoMaf['Tumor_Sample_Barcode'].isin(gliomaHyperIds)]
cohortMaf = mafAnnoMaf[mafAnnoMaf['Tumor_Sample_Barcode'].isin(colorectalHyperIds)]

#genes = set(['PTEN', 'PIK3CA', 'ARID1A', 'NF1', 'JAK1', 'RASA1', 'CIC', 'ATM'])
d = {}
cntr =0
for gene in genes:
    print cntr
    cntr += 1
    geneMaf = cohortMaf[cohortMaf['Hugo_Symbol'] == gene]
    l = []
    for case in set(geneMaf['Tumor_Sample_Barcode']):
        caseMafAllMuts = cohortMaf[cohortMaf['Tumor_Sample_Barcode'] == case]
        maxVaf = max(caseMafAllMuts['vaf'])
        caseGeneMaf = geneMaf[geneMaf['Tumor_Sample_Barcode'] == case]
        if caseGeneMaf.shape[0] > 1:
            vafs = list(caseGeneMaf['vaf'])
            minDist = min([abs(x - y) for i,x in enumerate(vafs) for j,y in enumerate(vafs) if i != j])
            minDistAsPct = minDist/maxVaf
            l.append(minDistAsPct)
    d[gene] = l
     
for gene in genes:
    if len(d[gene]) > 5:
        print gene, np.nanmedian(d[gene]), len(d[gene])   
        
###############################3333
print 


mutationMaf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/hypermutationAnalysisProjectData/hypermutationAnalysisProjectMaf.tsv')

print impactSigs[impactSigs['Nmut_Mb'] > 100].shape[0]/34744.0

impactSigs = mutationSigUtils.merge_signature_columns(impactSigs, mode='Stratton', drop=True, smokingMerge=True, confidence=True, mean=True, prefix='mean_')
impactSigs['dominantSignature'] = impactSigs.apply(lambda row: mutationSigUtils.get_dominant_signature(row.to_dict()), axis=1)

sigsToMakePlotWith = impactSigs[impactSigs['dominantSignature'].isin(set(['mean_1', 'mean_3', 'mean_SMOKING', 'mean_APOBEC', 'mean_7', 'mean_10', 'mean_MMR', 'mean_11', 'mean_14']))]
sigsToMakePlotWith = sigsToMakePlotWith[(sigsToMakePlotWith['Nmut'] > 10)]
sigsToMakePlotWith['orderingVal'] = sigsToMakePlotWith['dominantSignature'].apply(lambda x: 
    1 if x == 'mean_10'
    else 2 if x == 'mean_14'
    else 3 if x == 'mean_MMR'
    else 4 if x == 'mean_11'
    else 5 if x == 'mean_7'
    else 6 if x == 'mean_APOBEC'
    else 7 if x == 'mean_SMOKING'
    else 8
    )
sigsToMakePlotWith['dominantSignature'] = sigsToMakePlotWith['dominantSignature'].apply(lambda x:
    'POLE' if x == 'mean_10'
    else 'POLE & MMR' if x == 'mean_14'
    else 'MMR' if x == 'mean_MMR'
    else 'TMZ' if x == 'mean_11'
    else 'UV' if x == 'mean_7'
    else 'APOBEC' if x == 'mean_APOBEC'
    else 'SMOKING' if x == 'mean_SMOKING'
    else 'BRCA' if x == 'mean_3'
    else 'AGING' if x == 'mean_1'
    else x
    )
sigsToMakePlotWith.to_csv('~/Desktop/WORK/dataForLocalPlotting/signatureSummary.tsv', index=False, sep='\t')


################33
mmrDominant = set(impactSigs[(impactSigs['dominantSignature'] == 'mean_MMR') & (impactSigs['Nmut'] > 10)]['Tumor_Sample_Barcode'])
poleDominant = set(impactSigs[((impactSigs['dominantSignature'] == 'mean_10') | (impactSigs['dominantSignature'] == 'mean_14')) &(impactSigs['Nmut_Mb'] > 10)]['Tumor_Sample_Barcode'])
mutationMafMMR = mutationMaf[mutationMaf['Tumor_Sample_Barcode'].isin(mmrDominant)]
mutationMafMMR['Variant_Type']

dfBurdenInfo = maf_analysis_utils.summarize_per_case_mutation_info_for_mafs(mutationMaf)

dfBurdenInfo['oncogenicToNmutRatio'] = dfBurdenInfo.apply(lambda row: None if row['Nmut'] == 0 else 1.0*row['nOncogenicMuts']/row['Nmut'], axis=1)

l = []
for sigName in ['mean_MMR', 'mean_APOBEC', 'mean_10', 'mean_11', 'mean_7']:
    ids = set(impactSigs[((impactSigs['dominantSignature'] == sigName)) &(impactSigs['Nmut_Mb'] > 100)]['Tumor_Sample_Barcode'])
    ratio = np.nanmean(dfBurdenInfo[dfBurdenInfo['Tumor_Sample_Barcode'].isin(ids)]['oncogenicToNmutRatio'])
    caseMaf = mutationMaf[mutationMaf['Tumor_Sample_Barcode'].isin(ids)]
    indelFrac = 1.0*caseMaf[(caseMaf['Variant_Type'] == 'INS') | (caseMaf['Variant_Type'] == 'DEL')].shape[0]/caseMaf.shape[0]
    l.append({'sigName': sigName, 'ratio':ratio, 'indelFrac': indelFrac})
df = pd.DataFrame(l)
df.to_csv('~/Desktop/WORK/dataForLocalPlotting/indelAndSusceptibility.tsv', index=False, sep='\t')


print np.nanmean(dfBurdenInfo[dfBurdenInfo['Tumor_Sample_Barcode'].isin(poleDominant)]['oncogenicToNmutRatio'])



nMMRMutations = mutationMafMMR.shape[0]
nMMROncogenicMutations = mutationMafMMR[mutationMafMMR['oncogenic'].isin(oncogenicMutColNames)].shape[0]


print mutationMafMMR[(mutationMafMMR['Variant_Type'] == 'INS') | (mutationMafMMR['Variant_Type'] == 'DEL')].shape
print mutationMafMMR[(mutationMafMMR['Variant_Type'] == 'SNP')].shape




print mutationMafMMR[((mutationMafMMR['Variant_Type'] == 'INS') | (mutationMafMMR['Variant_Type'] == 'DEL')) & (mutationMafMMR['oncogenic'].isin(oncogenicMutColNames))].shape

print mutationMafMMR[((mutationMafMMR['Variant_Type'] == 'SNP')) & (mutationMafMMR['oncogenic'].isin(oncogenicMutColNames))].shape


print mutationMafMMR[mutationMafMMR['oncogenic'].isin(oncogenicMutColNames)]
