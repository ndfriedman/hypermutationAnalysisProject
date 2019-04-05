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
import variant_consequence_analysis_util

def find_pan_impact_mut_susceptibility_prob(df):
    refTrisWhereStopCodonCanOccur = variant_consequence_analysis_util.get_all_ref_tris_where_stop_codon_can_occur()
    poleMotifMuts = ['TCT', 'TCG']
    return variant_consequence_analysis_util.analyze_codon_prevalence(df, poleMotifMuts)

def summarize_gene_stop_codon_content(df, analyzeMaf, motifs):
    listOfDicts = []
    tumorSupressors = analysis_utils.get_tumor_supressor_genes()
    for gene in tumorSupressors:
        geneSeqDf = df[df['gene'] == gene]
        n = find_pan_impact_mut_susceptibility_prob(geneSeqDf)
        truncatingDf = analyzeMaf[(analyzeMaf['Hugo_Symbol'] == gene) & (analyzeMaf['Consequence'] == 'stop_gained')]
        truncatingDfSpecificMotifs = truncatingDf[truncatingDf['Ref_Tri'].isin(motifs)]
        nObservedTruncatingInGene = truncatingDf.shape[0]
        nTruncatingDfSpecificMotifs = truncatingDfSpecificMotifs.shape[0]
        listOfDicts.append({'Gene': gene, 'NPoleStopCodonPossible': n, 'NObservedTrunc': nObservedTruncatingInGene, 'NObsercedTruncAtMotif': nTruncatingDfSpecificMotifs})
    return pd.DataFrame(listOfDicts)

sequenceInfoDf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/impactGeneSequenceSummary.tsv')

maf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/hypermutationAnalysisProjectData/hypermutationAnalysisProjectMaf.tsv')
sigs = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/hypermutationAnalysisProjectData/hypermutationAnalysisProjectSignatures.tsv')

sigs['pid'] = sigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
sigs['cancer_type'] = sigs['pid'].apply(lambda x: cDict[x] if x in cDict else None)
maf['quadNuc'] = maf.apply(lambda row: mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)
    
hypermuationThreshold = 80
hypermutatorIds = set(sigs[sigs['Nmut_Mb'] > hypermuationThreshold]['Tumor_Sample_Barcode'])
maf['isHypermutator'] = maf['Tumor_Sample_Barcode'].apply(lambda x: True if x in hypermutatorIds else False)
endometrialIds = set(sigs[sigs['cancer_type'] == 'Endometrial Cancer']['Tumor_Sample_Barcode'])
endometrialHypermutatorMaf = maf[(maf['Tumor_Sample_Barcode'].isin(endometrialIds)) & (maf['isHypermutator'] == True)]
endometrialHypermutatorMafOncogenic = endometrialHypermutatorMaf[endometrialHypermutatorMaf['oncogenic'].notnull()]
endometrialHypermutatorMafOncogenic['quadNuc'] = endometrialHypermutatorMafOncogenic.apply(lambda row: mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)
poleSigIds = set(sigs[(sigs['mean_10'] > .15)]['Tumor_Sample_Barcode'])
endometrialPoleMaf = endometrialHypermutatorMaf[endometrialHypermutatorMaf['Tumor_Sample_Barcode'].isin(poleSigIds)]

endometrialPoleMaf[endometrialPoleMaf['quadNuc']]

print 1.0*sigs[sigs['Nmut_Mb'] > 80].shape[0]/sigs.shape[0]

#ptenOnlyDf = sequenceInfoDf[sequenceInfoDf['gene'] == 'PTEN']
#find_pan_impact_mut_susceptibility_prob(ptenOnlyDf)

print maf[(maf['Tumor_Sample_Barcode'] == 'P-0028776-T01-IM6')&(maf['Hugo_Symbol'] == 'LATS2')]['HGVSp']

allCodons = variant_consequence_analysis_util.get_all_ref_tris_where_stop_codon_can_occur()
df = summarize_gene_stop_codon_content(sequenceInfoDf, endometrialPoleMaf, allCodons)
df.to_csv('~/Desktop/WORK/dataForLocalPlotting/truncatingInfoPole.tsv', index=False, sep='\t')
#################################3

case = 'P-0028776-T01-IM6'
gene = 'MSH2'
positions = list(endometrialHypermutatorMaf[(endometrialHypermutatorMaf['Tumor_Sample_Barcode'] == case) & (endometrialHypermutatorMaf['Hugo_Symbol'] == gene)]['Start_Position'])
for position in positions:
    for p in positions:
        print abs(position - p)



print endometrialHypermutatorMaf[(endometrialHypermutatorMaf['Tumor_Sample_Barcode'] == case) & (endometrialHypermutatorMaf['Hugo_Symbol'] == gene)]['Start_Position']
print endometrialHypermutatorMaf[(endometrialHypermutatorMaf['Tumor_Sample_Barcode'] == case) & (endometrialHypermutatorMaf['Hugo_Symbol'] == gene)]['HGVSp']




Counter(df['Mutation Effect'])
df['Protein Change']
Counter(df['Oncogenicity'])

impactGenes = set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KMT2B', 'KMT5A', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51C', 'RAD51L1', 'RAD51L3', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2'])

endometrialMafAnnoMaf['vaf'] = endometrialMafAnnoMaf.apply(lambda row: 0 if row['t_ref_count'] == 0 else 1.0*row['t_alt_count']/row['t_ref_count'], axis=1)

endometrialMafAnnoMaf[(endometrialMafAnnoMaf['Tumor_Sample_Barcode'] == 'P-0011569-T01-IM5') & (endometrialMafAnnoMaf['ccf_1copy_lower'] >.99)]['Hugo_Symbol']

for case in set(endometrialMafAnnoMaf['Tumor_Sample_Barcode']):
    caseMaf = endometrialMafAnnoMaf[endometrialMafAnnoMaf['Tumor_Sample_Barcode'] == case]
    for gene in impactGenes:
        geneMaf = caseMaf[caseMaf['Hugo_Symbol'] == gene]
        snps = geneMaf[geneMaf['Variant_Type'] == 'SNP']
        positions = geneMaf['Start_Position']
        positionsCopy = positions
        for position in positions:
            for p in positionsCopy:
                diff = abs(position - p)
                if diff < 100 and diff != 0:
                    print case, gene
                    print geneMaf['ccf_1copy_lower']
                    #if geneMaf[geneMaf['ccf_1copy_lower'] > .99].shape[0] >= 1:
                    #    print gene, diff, geneMaf['ccf_1copy_lower']
                    



######################################################################33

mafAnnoMaf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/All.dmp_somatic_data_mutations_unfiltered.mafAnno.oncokb.hotspots.maf', sep='\t')
mafAnnoMaf['isHypermutator'] = mafAnnoMaf['Tumor_Sample_Barcode'].apply(lambda x: True if x in hypermutatorIds else False)

endometrialMafAnnoMaf = mafAnnoMaf[mafAnnoMaf['Tumor_Sample_Barcode'].isin(endometrialIds)]
endometrialHypermuatorMafAnnoMaf = endometrialMafAnnoMaf[endometrialMafAnnoMaf['isHypermutator'] == True]
endometrialNotHypermuatorMafAnnoMaf = endometrialMafAnnoMaf[endometrialMafAnnoMaf['isHypermutator'] == False]
endometrialHypermuatorMafAnnoMafPten = endometrialHypermuatorMafAnnoMaf[endometrialHypermuatorMafAnnoMaf['Hugo_Symbol'] == 'PTEN']
endometrialNotHypermuatorMafAnnoMafPten = endometrialNotHypermuatorMafAnnoMaf[endometrialNotHypermuatorMafAnnoMaf['Hugo_Symbol'] == 'PTEN']


endometrialHypermuatorMafAnnoMafPten[endometrialHypermuatorMafAnnoMafPten['lcn'] == 1]['Tumor_Sample_Barcode']

len(set(endometrialNotHypermuatorMafAnnoMafPten[endometrialNotHypermuatorMafAnnoMafPten['lcn'] == 1]['Tumor_Sample_Barcode']))

for case in set(endometrialNotHypermuatorMafAnnoMafPten['Tumor_Sample_Barcode']):
    caseMaf = endometrialNotHypermuatorMafAnnoMafPten[endometrialNotHypermuatorMafAnnoMafPten['Tumor_Sample_Barcode'] == case]
    if caseMaf.shape[0] == 1:
        print caseMaf['lcn'].iloc[0]
    

###############################################################################
df = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/allAnnotatedVariantsOncokb.txt')
df.columns.values

df['Protein Change']

print len(set(maf['Tumor_Sample_Barcode']))

print len(set(maf[maf['Tumor_Sample_Barcode'] == '']))








    
    



