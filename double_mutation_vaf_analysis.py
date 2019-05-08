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
import clonality_analysis_util

#todo build data across cohort case by case
#note there are different analyses including matrix of relevant/relevant and gene specific bialleic 
#compare vafs between pairs of different classes of mutations in a single case
def compare_vafs(maf, referenceCohortMaf):
    
    returnDict = {'allMutations': [], 'allOncogenicMutations':[],
                  'allRecurrentTumorSuppressors': [], 'oncogenicRecurrentTumorSuppresors': [],
                  'allRecurrentOncogenes': [], 'oncogenicRecurrentOncogenes': []}
    
    recurrentTumorSupressors, recurrentOncogenes = maf_analysis_utils.enumerate_recurrently_mutated_tumor_supressors_and_oncogenes(referenceCohortMaf, thresh=.05)
    oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
    cases = set(maf['Tumor_Sample_Barcode'])
    genes = set(maf['Hugo_Symbol'])
    
    print recurrentOncogenes
    
    for case in cases:
        print case
        #for gene in genes:
         #   geneMaf = maf[maf['Hugo_Symbol'] == gene]
            #here are all the cases/areas we look at
        caseMafAllMuts = maf[maf['Tumor_Sample_Barcode'] == case]
        #caseMafAllMuts = geneMaf[geneMaf['Tumor_Sample_Barcode'] == case]
        
        caseMafAllOncogenicMuts = caseMafAllMuts[caseMafAllMuts['oncogenic'].isin(oncogenicMutColNames)]
        caseMafAllReccurentTumorSuppressors = caseMafAllMuts[caseMafAllMuts['Hugo_Symbol'].isin(recurrentTumorSupressors)]
        caseMafOncogenicRecurrentTumorSuppresors = caseMafAllReccurentTumorSuppressors[caseMafAllReccurentTumorSuppressors['oncogenic'].isin(oncogenicMutColNames)]
        caseMafAllRecurrentOncogenes = caseMafAllMuts[caseMafAllMuts['Hugo_Symbol'].isin(recurrentOncogenes)]
        caseMafOncogenicRecurrentOncogenes = caseMafAllRecurrentOncogenes[caseMafAllRecurrentOncogenes['oncogenic'].isin(oncogenicMutColNames)]
        
        #add information to our return dict
        returnDict['allMutations'] = returnDict['allMutations'] + analysis_utils.calculate_all_pairwise_differences(np.array(list(caseMafAllMuts['ccf_Mcopies'])))
        returnDict['allOncogenicMutations'] = returnDict['allOncogenicMutations'] + analysis_utils.calculate_all_pairwise_differences(np.array(list(caseMafAllOncogenicMuts['ccf_Mcopies'])))
        returnDict['allRecurrentTumorSuppressors'] = returnDict['allRecurrentTumorSuppressors'] + analysis_utils.calculate_all_pairwise_differences(np.array(list(caseMafAllReccurentTumorSuppressors['ccf_Mcopies'])))
        returnDict['oncogenicRecurrentTumorSuppresors'] = returnDict['oncogenicRecurrentTumorSuppresors'] + analysis_utils.calculate_all_pairwise_differences(np.array(list(caseMafOncogenicRecurrentTumorSuppresors['ccf_Mcopies'])))
        returnDict['allRecurrentOncogenes'] = returnDict['allRecurrentOncogenes'] + analysis_utils.calculate_all_pairwise_differences(np.array(list(caseMafAllRecurrentOncogenes['ccf_Mcopies'])))
        returnDict['oncogenicRecurrentOncogenes'] = returnDict['oncogenicRecurrentOncogenes'] + analysis_utils.calculate_all_pairwise_differences(np.array(list(caseMafOncogenicRecurrentOncogenes['ccf_Mcopies'])))
        
    return returnDict
       
        
mafAnnoMaf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/All.dmp_somatic_data_mutations_unfiltered.mafAnno.oncokb.hotspots.maf', sep='\t')
mafAnnoMaf['vaf'] = mafAnnoMaf.apply(lambda row: 1.0*row['t_alt_count']/(row['t_alt_count'] + row['t_ref_count']) if row['t_alt_count'] + row['t_ref_count'] > 0 else None, axis = 1)

impactSigs = pd.read_table(pathPrefix + '/ifs/res/taylorlab/impact_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
impactSigs['pid'] = impactSigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
impactSigs['cancer_type'] = impactSigs['pid'].apply(lambda x: cDict[x] if x in cDict else None)

hypermutationThresh = 50
endometrialHyperIds = set(impactSigs[(impactSigs['cancer_type'] == 'Endometrial Cancer') & (impactSigs['Nmut_Mb'] > hypermutationThresh)]['Tumor_Sample_Barcode'])
endometrialNotHyperIds = set(impactSigs[(impactSigs['cancer_type'] == 'Endometrial Cancer') & (impactSigs['Nmut_Mb'] <= hypermutationThresh)]['Tumor_Sample_Barcode'])
gliomaHyperIds = set(impactSigs[(impactSigs['cancer_type'] == 'Glioma') & (impactSigs['Nmut_Mb'] > hypermutationThresh)]['Tumor_Sample_Barcode'])
colorectalHyperIds = set(impactSigs[(impactSigs['cancer_type'] == 'Colorectal Cancer') & (impactSigs['Nmut_Mb'] > hypermutationThresh)]['Tumor_Sample_Barcode'])
colorectalNotHyperIds = set(impactSigs[(impactSigs['cancer_type'] == 'Colorectal Cancer') & (impactSigs['Nmut_Mb'] <= hypermutationThresh)]['Tumor_Sample_Barcode'])

genes = set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KMT2B', 'KMT5A', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51C', 'RAD51L1', 'RAD51L3', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2'])
oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])

mafAnnoMafWithData = mafAnnoMaf[mafAnnoMaf['ccf_Mcopies'].notnull()]
mafAnnoMafWithData['pid'] = mafAnnoMafWithData['Tumor_Sample_Barcode'].apply(lambda x: x[:9])

#cohortMaf = mafAnnoMaf[mafAnnoMaf['Tumor_Sample_Barcode'].isin(endometrialHyperIds) & (mafAnnoMaf['oncogenic'].isin(oncogenicMutColNames))]
#cohortMaf = mafAnnoMaf[mafAnnoMaf['Tumor_Sample_Barcode'].isin(endometrialHyperIds)]
#cohortMaf = mafAnnoMaf[mafAnnoMaf['Tumor_Sample_Barcode'].isin(gliomaHyperIds)]
cohortMaf = mafAnnoMafWithData[mafAnnoMafWithData['Tumor_Sample_Barcode'].isin(endometrialHyperIds)]


#####
#########
##############
#################   NEW VERSION   ########### 

#TODO check for genome doubling
def create_two_hit_mutation_summary(df, alleleCountThresh =5, doubleThresh=1.75):
    
    listOfDicts = []
    #PROCESS THE MAF AS NEEDED
    df = maf_analysis_utils.fix_mll_genes(df)
    df = maf_analysis_utils.mark_cases_with_flat_genomes(df)
    df = maf_analysis_utils.mark_cases_with_median_vaf_of_case(df)
    df = df[df['oncogenic'].notnull()]
    df = df[df['HGVSp_Short'].notnull()]
    
    df['isDouble'] = df.apply(lambda row: clonality_analysis_util.is_mut_double_hit(row, row['medianVaf'], row['FlatGenome'], doubleThresh), axis=1)
    
    dictOfMuts = dict(Counter(df[df['isDouble'] == True]['Hugo_Symbol']))
    for gene, count in dictOfMuts.items():
        if count > 1:
            for allele, aCount in Counter(df[(df['Hugo_Symbol'] == gene) & (df['isDouble'] == True)]['HGVSp_Short']).items():
                indel = False
                if '*' in allele: indel = True
                d = {'Gene': gene, 'Allele': allele, 'AlleleCount': aCount, 'count': count, 'orderingVal': count, 'indel':indel}
                listOfDicts.append(d)
        else: #PROCESS DATA SEPARATELY FOR ONE OFF ALELLES
            allele = df[(df['Hugo_Symbol'] == gene) & (df['isDouble'] == True)].iloc[0]['HGVSp_Short']
            indel = False
            if '*' in allele: indel = True
            d = {'Gene': 'other', 'Allele': allele, 'AlleleCount': 1, 'count': 1, 'orderingVal': -1, 'indel':indel}
            listOfDicts.append(d)
    
    df = pd.DataFrame(listOfDicts)
    df['label'] = df.apply(lambda row: row['Allele'] if row['AlleleCount'] >= alleleCountThresh and row['Allele'] != 'various' else None, axis=1)
    return df
    
#ANALYZE COHORTS TO SEE PREVALENCE AND SPECIFICITY OF DOUBLE HIT MUTS

colorectalDf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/subsettedMafs/Colorectal_HypermutantCaseMuts_MAF_ANNO.maf')
endometrialDf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/subsettedMafs/Endometrial_HypermutantCaseMuts_MAF_ANNO.maf')
gliomaDf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/subsettedMafs/Glioma_HypermutantCaseMuts_MAF_ANNO.maf')

#TCF7L2--9mer, RNF43---7mer MSH3 AAAAAAA homopolymer, PTEN 6 homopolymer
    #APC  p.T1556Nfs*3--6mer, p.R856Nfs*6: 3mer

endoSummary = create_two_hit_mutation_summary(endometrialDf)
endoSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/endometrialDoubleHitPlotting.tsv', index=False, sep='\t')

colorectalSummary =  create_two_hit_mutation_summary(colorectalDf)
colorectalSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/colorectalDoubleHitPlotting.tsv', index=False, sep='\t')

gliomaSummary =  create_two_hit_mutation_summary(gliomaDf)
gliomaSummary.to_csv('~/Desktop/WORK/dataForLocalPlotting/gliomaDoubleHitPlotting.tsv', index=False, sep='\t')

print 

################3
gliomaDf2 = maf_analysis_utils.mark_cases_with_median_vaf_of_case(gliomaDf)
gliomaDf2 = gliomaDf2[gliomaDf2['oncogenic'].notnull()]
gliomaDf2 = gliomaDf2[gliomaDf2['HGVSp_Short'].notnull()]
gliomaDf2['isDouble'] = gliomaDf2.apply(lambda row: clonality_analysis_util.is_mut_double_hit(row, row['medianVaf'], row['FlatGenome'], 1.75), axis=1)
    
for i in gliomaDf2[(gliomaDf2['isDouble'] == True) & (gliomaDf['Hugo_Symbol'] == 'IDH1')]['t_var_freq']:
    print i


