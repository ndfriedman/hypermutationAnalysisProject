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

#returns a dataframe mapping case names, n oncogenic mutations and n hotspots, fraction hotspots at signature enriched motif per case
def enumerate_case_mutation_info_summary(df, enrichedSigMotifs):
    listOfDicts = []
    cases = set(df['Tumor_Sample_Barcode'])
    cntr = 0
    for case in cases:
        if cntr%100 == 0: 
            print cntr, len(cases)
        cntr += 1
        
        localD = {}
        caseDf = df[df['Tumor_Sample_Barcode'] == case]
        
        #TODO EXTEND THIS ANALYSIS TO INCORP FRAC HOTSPOTS at enriched motif and FRAC oncogenic muts at enriched MOTIF
        
        oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic']) #enumerate col names for likely oncogenic mutations
        
        hotspotMutDf = caseDf[caseDf['is-a-hotspot'] == 'Y']
        oncogenicMutDf = caseDf[caseDf['oncogenic'].isin(oncogenicMutColNames)]
        oncogenicOrHotspotMutations = caseDf[(caseDf['oncogenic'].isin(oncogenicMutColNames)) | (caseDf['is-a-hotspot'] == 'Y')]
        
        nHotspotMutations = hotspotMutDf.shape[0]
        nOncogenicMutations = oncogenicMutDf.shape[0]
        nOncogenicOrHotspotMutations = oncogenicOrHotspotMutations.shape[0] #we need to count this separately because they may overlap 
        
        fracOncogenicMutationsAtEnrichedMotif = None
        fracHotpsotMutationsAtEnrichedMotif = None
        fracDriverMutationsAtEnrichedMotif = None
        if nHotspotMutations > 0:
            fracHotpsotMutationsAtEnrichedMotif = 1.0*hotspotMutDf[hotspotMutDf['quadNuc'].isin(enrichedSigMotifs)].shape[0]/nHotspotMutations
        if nOncogenicOrHotspotMutations > 0:
            fracOncogenicMutationsAtEnrichedMotif = 1.0*oncogenicMutDf[oncogenicMutDf['quadNuc'].isin(enrichedSigMotifs)].shape[0]/nOncogenicMutations
        if nOncogenicOrHotspotMutations > 0:
            fracDriverMutationsAtEnrichedMotif = 1.0*oncogenicOrHotspotMutations[oncogenicOrHotspotMutations['quadNuc'].isin(enrichedSigMotifs)].shape[0]/nOncogenicOrHotspotMutations
        
        nMut = caseDf.shape[0]
        nMutToHotspotRatio = None
        nMutToOncogenicRatio = None
        nMutToOncogenicAndHotspotRatio = None
        if nHotspotMutations > 0:
            nMutToHotspotRatio = 1.0*nMut/nHotspotMutations
        if nOncogenicMutations > 0:
            nMutToOncogenicRatio = 1.0*nMut/nOncogenicMutations
        if nOncogenicOrHotspotMutations > 0:
            nMutToOncogenicAndHotspotRatio = 1.0*nMut/nOncogenicOrHotspotMutations
        
        
        #add in all the information to the local dict
        localD['Tumor_Sample_Barcode'] = case
        localD['nHotspots'] = nHotspotMutations
        localD['nOncogenicMutations'] = nOncogenicMutations
        localD['nOncogenicOrHotspotMutations'] = nOncogenicOrHotspotMutations
        localD['fracOncogenicMutationsAtEnrichedMotif'] = fracOncogenicMutationsAtEnrichedMotif
        localD['fracHotspotMutationsAtEnrichedMotif'] = fracHotpsotMutationsAtEnrichedMotif
        localD['fracDriverMutationsAtEnrichedMotif'] = fracDriverMutationsAtEnrichedMotif
        localD['nMutToHotspotRatio'] = nMutToHotspotRatio
        localD['nMutToOncogenicRatio'] = nMutToOncogenicRatio
        localD['nMutToOncogenicAndHotspotRatio'] = nMutToOncogenicAndHotspotRatio
        localD['Nmut'] = nMut
        
        listOfDicts.append(localD)
      
    df = pd.DataFrame(listOfDicts)
    return df

mafWithInfo = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/annotatedOncoPlusHotspotMafAllImpact_trinuc')
impactSigs = pd.read_table(pathPrefix + '/ifs/res/taylorlab/impact_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
impactSigs['pid'] = impactSigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
impactSigs['cancer_type'] = impactSigs['pid'].apply(lambda x: cDict[x] if x in cDict else None)
mafWithInfo['pid'] = mafWithInfo['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
mafWithInfo['cancer_type'] = mafWithInfo['pid'].apply(lambda x: cDict[x] if x in cDict else None)

#adjust column names to make the 
renameDict = {key:value for (key,value) in [('mean_' + str(i), 'Signature.' + str(i)) for i in range(1,31)]}
impactSigs = impactSigs.rename(columns=renameDict)
impactSigs = mutationSigUtils.merge_signature_columns(impactSigs, mode='Stratton', drop=False, smokingMerge=True, confidence=False, mean=True, prefix='Signature.')

nmut_MbDict = dict(zip(impactSigs['Tumor_Sample_Barcode'], impactSigs['Nmut_Mb']))

cancerTypesToFocusOn = set(['Colorectal Cancer', 'Non-Small Cell Lung Cancer', 'Glioma', 'Melanoma', 'Endometrial Cancer', 'Bladder Cancer'])

#mafLimitedToCancerType = mafWithInfo[mafWithInfo['cancer_type'].isin(cancerTypesToFocusOn)]

#ADD IN INFORMATION ABOUT QUADNUCs
#mafLimitedToCancerType['quadNuc'] = mafLimitedToCancerType.apply(lambda row: mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)
mafWithInfo['quadNuc'] = mafWithInfo.apply(lambda row: mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)

#infoDf  = enumerate_case_mutation_info_summary(mafLimitedToCancerType, set([]))
infoDf  = enumerate_case_mutation_info_summary(mafWithInfo, set([]))

infoDf['Nmut_Mb'] = infoDf['Tumor_Sample_Barcode'].apply(lambda x:  nmut_MbDict[x] if x in nmut_MbDict else None)
infoDf['pid'] = infoDf['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
infoDf['cancer_type'] = infoDf['pid'].apply(lambda x: cDict[x] if x in cDict else None)

#CREATE THE COHORTS!
signatureThreshold = .2 #the threshold at which we call a signature as existing
notHighThresh = 15
highThresh = 50

#POLE is divided into two classes: POLE endometrial and other pole
#If a case is 'mixed' POLE plus MMR I consider it POLE
poleEndometrialIds = set(impactSigs[(impactSigs['cancer_type'] == 'Endometrial Cancer') & (impactSigs['Nmut_Mb'] > highThresh) 
& ((impactSigs['Signature.10'] > signatureThreshold) | (impactSigs['Signature.14'] > signatureThreshold))]['Tumor_Sample_Barcode'])
    
otherPoleIds = set(impactSigs[(impactSigs['cancer_type'] != 'Endometrial Cancer') & (impactSigs['Nmut_Mb'] > highThresh) 
& ((impactSigs['Signature.10'] > signatureThreshold) | (impactSigs['Signature.14'] > signatureThreshold))]['Tumor_Sample_Barcode'])

#MMR is divided into three classes: Colorectal MMR, Endometrial MMR and Other MMR
colorectalMMRIds = set(impactSigs[(impactSigs['cancer_type'] == 'Colorectal Cancer') & (impactSigs['Signature.MMR'] > signatureThreshold) 
& ((impactSigs['Signature.10'] < signatureThreshold) & (impactSigs['Signature.14'] < signatureThreshold))]['Tumor_Sample_Barcode'])

endometrialMMRIds = set(impactSigs[(impactSigs['cancer_type'] == 'Endometrial Cancer') & (impactSigs['Signature.MMR'] > signatureThreshold) 
& ((impactSigs['Signature.10'] < signatureThreshold) & (impactSigs['Signature.14'] < signatureThreshold))]['Tumor_Sample_Barcode'])
 
otherMMRIds = set(impactSigs[(impactSigs['cancer_type'] != 'Colorectal Cancer') & (impactSigs['cancer_type'] != 'Endometrial Cancer') & (impactSigs['Nmut_Mb'] > highThresh)  & (impactSigs['Signature.MMR'] > signatureThreshold)
& ((impactSigs['Signature.10'] < signatureThreshold) | (impactSigs['Signature.14'] < signatureThreshold))]['Tumor_Sample_Barcode'])

#TMZ hypermutated is a specific group
gliomaTMZIds = set(impactSigs[(impactSigs['cancer_type'] == 'Glioma') & (impactSigs['Nmut_Mb'] > highThresh) 
& (impactSigs['Signature.11'] > signatureThreshold)]['Tumor_Sample_Barcode'])

#other High groups are just based on TMB and cancer type
cancerTypesForControls = set(['Colorectal Cancer', 'Non-Small Cell Lung Cancer', 'Glioma', 'Melanoma', 'Endometrial Cancer', 'Bladder Cancer'])    

infoDf['cohort'] = infoDf.apply(lambda row:
    
    row['cancer_type'] + '__not_high' if row['Nmut_Mb'] < notHighThresh and row['cancer_type'] in cancerTypesForControls
    else 'POLE_Endometrial' if row['Tumor_Sample_Barcode'] in poleEndometrialIds
    else 'POLE_Other' if row['Tumor_Sample_Barcode'] in otherPoleIds
    else 'MMR_Other' if row['Tumor_Sample_Barcode'] in otherMMRIds
    else 'MMR_Colorectal' if row['Tumor_Sample_Barcode'] in colorectalMMRIds
    else 'MMR_Endometrial' if row['Tumor_Sample_Barcode'] in endometrialMMRIds
    else 'MMR_Other' if row['Tumor_Sample_Barcode'] in otherMMRIds
    else 'TMZ_glioma' if row['Tumor_Sample_Barcode'] in gliomaTMZIds
    else row['cancer_type'] + '_high' if row['Nmut_Mb'] > highThresh and row['cancer_type'] in set(['Non-Small Cell Lung Cancer', 'Melanoma', 'Bladder Cancer'])
    else None
    ,axis=1)
 
#OLD WAY    
#highMutBurdenThresh = 20
#hypermutatorMutBurdenThresh = 80
#infoDf['cohort'] = infoDf.apply(lambda row: 
#    row['cancer_type'] + '__not_high' if row['Nmut_Mb'] < 10
#    else row['cancer_type'] + '__high' if row['Nmut_Mb'] < 80
#    else row['cancer_type'] + '__hypermutant' if row['Nmut_Mb'] >= 80
#    else row['cancer_type'] + '__not_high' #simplify because almost all cases without an nmut_mb estimate are not hypermutated #TODO FIX
#    , axis = 1)
#orderingDict1 = {'Endometrial Cancer': 0, 'Colorectal Cancer': 1, 'Glioma': 2, 'Melanoma': 3, 'Bladder Cancer': 4, 'Non-Small Cell Lung Cancer':5}
#orderingDict2 = {'hypermutant': .1, 'high': .2, 'not_high': .3}
#infoDf['orderingVal'] = infoDf['cohort'].apply(lambda x: orderingDict1[x.split('__')[0]] + orderingDict2[x.split('__')[1]]) 

orderingValDict = {'POLE_Other': 0, 'MMR_Other': 1,
                   'POLE_Endometrial':2, 'MMR_Endometrial':3, 'Endometrial Cancer__not_high':4,
                   'MMR_Colorectal':5, 'Colorectal Cancer__not_high':6,
                   'TMZ_glioma':7, 'Glioma__not_high':8,
                   'Melanoma_high':9, 'Melanoma__not_high':10,
                   'Bladder Cancer_high':11, 'Bladder Cancer__not_high':12,
                   'Non-Small Cell Lung Cancer_high':13, 'Non-Small Cell Lung Cancer__not_high':14,
                   }

infoDf['orderingVal'] = infoDf['cohort'].apply(lambda x: orderingValDict[x] if x in orderingValDict else None)
    

#NO LONGER NEED TO DO THIS COHORTS HAVE THIS BAKED IN
###ADD DOMINANT SIGNATURE INFORMATION
#sigNames = [i for i in impactSigs.columns.values if 'Signature.' in i]
#impactSigs['dominantSignature'] = impactSigs.apply(lambda row: mutationSigUtils.get_dominant_signature(row.to_dict(), cols=sigNames, prefix='Signature.'), axis=1)
#dominantSigDict = dict(zip(impactSigs['Tumor_Sample_Barcode'], impactSigs['dominantSignature']))
#infoDf['dominantSignature'] = infoDf['Tumor_Sample_Barcode'].apply(lambda x: dominantSigDict[x] if x in dominantSigDict else None)
#infoDf['dominantSignature'] = infoDf.apply(lambda row: 'Not Enough Mutations' if row['Nmut'] < 10 else row['dominantSignature'], axis=1)  
#infoDf['signatureAetiology'] = infoDf['dominantSignature'].apply(lambda x:
#    'Smoking' if x == 'Signature.SMOKING'
#    else 'MMR' if x == 'Signature.MMR'
#    else 'APOBEC' if x == 'Signature.APOBEC'
#    else 'Age' if x == 'Signature.1'
#    else 'POLE' if x == 'Signature.10'
#    else 'Mixed POLE/MMR' if x == 'Signature.14'
#    else 'UV' if x == 'Signature.7'
#    else 'BRCA' if x == 'Signature.3'
#    else 'TMZ' if x == 'Signature.11'
#    else x if x == 'Not Enough Mutations'
#    else 'Other Signature'
#    )

#Age of diagnosis stuff
ageAtSequencing = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/clinical_msk_impact_repository_data_adj.txt')
ageAtSequencingDict = dict(zip(ageAtSequencing['#Sample Identifier'], ageAtSequencing['Age at Which Sequencing was Reported (Days)']))
    
infoDf['ageAtSequencing'] = infoDf['Tumor_Sample_Barcode'].apply(lambda x: ageAtSequencingDict[x] if x in ageAtSequencingDict else None)

###ADD more information about doublets
#summaryDfDoubleMuts = double_mutation_analysis_util.create_double_mutation_summary_vanilla_maf(mafLimitedToCancerType)
#nGeneDoubleWithOncogenicDict = dict(zip(summaryDfDoubleMuts['Tumor_Sample_Barcode'], summaryDfDoubleMuts['nGenesDoubleOncogenicPerCase']))
#infoDf['nGenesWithDoubleOncogenic'] = infoDf['Tumor_Sample_Barcode'].apply(lambda x: nGeneDoubleWithOncogenicDict[x] if x in nGeneDoubleWithOncogenicDict else None) 

infoDf['cancer_type_fill'] = infoDf.apply(lambda row: 'Other' if row['cohort'] == 'MMR_Other' or row['cohort'] == 'POLE_Other' else row['cancer_type'], axis=1)
infoDf[infoDf['cohort'].notnull()].to_csv('~/Desktop/WORK/dataForLocalPlotting/mutburdenBoxplotV2.tsv', sep='\t', index=False)



##################SIGNATURE COHORT INFO


impactSigs = pd.read_table(pathPrefix + '/ifs/res/taylorlab/impact_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
impactSigs['pid'] = impactSigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
renameDict = {key:value for (key,value) in [('mean_' + str(i), 'Signature.' + str(i)) for i in range(1,31)]}
impactSigs = impactSigs.rename(columns=renameDict)
impactSigs = mutationSigUtils.merge_signature_columns(impactSigs, mode='Stratton', drop=False, smokingMerge=True, confidence=False, mean=True, prefix='Signature.')

sigNames = [i for i in impactSigs.columns.values if 'Signature.' in i]
impactSigs['dominantSignature'] = impactSigs.apply(lambda row: mutationSigUtils.get_dominant_signature(row.to_dict(), cols=sigNames, prefix='Signature.'), axis=1)
impactSigs['cancer_type'] = impactSigs['pid'].apply(lambda x: cDict[x] if x in cDict else None)
dominantSigDict = dict(zip(impactSigs['Tumor_Sample_Barcode'], impactSigs['dominantSignature']))
impactSigs['dominantSignature'] = impactSigs.apply(lambda row: 'Not Enough Mutations' if row['Nmut'] < 10 else row['dominantSignature'], axis=1)  
impactSigs['signatureAetiology'] = impactSigs['dominantSignature'].apply(lambda x:
    'Smoking' if x == 'Signature.SMOKING'
    else 'MMR' if x == 'Signature.MMR'
    else 'APOBEC' if x == 'Signature.APOBEC'
    else 'Age' if x == 'Signature.1'
    else 'POLE' if x == 'Signature.10'
    else 'Mixed POLE/MMR' if x == 'Signature.14'
    else 'UV' if x == 'Signature.7'
    else 'BRCA' if x == 'Signature.3'
    else 'TMZ' if x == 'Signature.11'
    else x if x == 'Not Enough Mutations'
    else 'Other Signature'
    )
    
impactSigs.to_csv('~/Desktop/WORK/dataForLocalPlotting/signaturePlotting.tsv', index=False, sep='\t')



###
###
#########
##############
####################

#RANDOM WORKSPACE

endometrialHypermutators = set(impactSigs[(impactSigs['cancer_type'] == 'Endometrial Cancer') & (impactSigs['Nmut_Mb'] > 80)]['Tumor_Sample_Barcode'])
eHyperMaf = mafWithInfo[mafWithInfo['Tumor_Sample_Barcode'].isin(endometrialHypermutators)]

colorectalHypermutators = set(impactSigs[(impactSigs['cancer_type'] == 'Colorectal Cancer') & (impactSigs['Nmut_Mb'] > 50)]['Tumor_Sample_Barcode'])
cHyperMaf = mafWithInfo[mafWithInfo['Tumor_Sample_Barcode'].isin(colorectalHypermutators)]

gliomaHypermutators = set(impactSigs[(impactSigs['cancer_type'] == 'Glioma') & (impactSigs['Nmut_Mb'] > 50)]['Tumor_Sample_Barcode'])
gHyperMaf = mafWithInfo[mafWithInfo['Tumor_Sample_Barcode'].isin(gliomaHypermutators)]

bladderHypermutators = set(impactSigs[(impactSigs['cancer_type'] == 'Bladder Cancer') & (impactSigs['Nmut_Mb'] > 50)]['Tumor_Sample_Barcode'])
bHyperMaf = mafWithInfo[mafWithInfo['Tumor_Sample_Barcode'].isin(bladderHypermutators)]


analyzeMaf = cHyperMaf

oncogenicMutCols = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
impact468 = set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'AMER1', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'GTF2I', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MLL4', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SETD8', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2']) 
listOfDicts = []
for case in set(analyzeMaf['Tumor_Sample_Barcode']):
    print case
    localD = {}
    #localD = {'Tumor_Sample_Barcode': case}
    for gene in impact468:
        caseMaf = analyzeMaf[analyzeMaf['Tumor_Sample_Barcode'] == case]
        nGeneMuts = caseMaf[(caseMaf['Hugo_Symbol'] == gene) & caseMaf['oncogenic'].isin(oncogenicMutCols)].shape[0]
        localD['n_' + gene + '_muts'] = nGeneMuts
    print Counter(localD.values())
    listOfDicts.append(localD)
df = pd.DataFrame(listOfDicts)
df.to_csv('~/Desktop/WORK/dataForLocalPlotting/mutationNumberDistributions_bladder.tsv', index=False, sep='\t')

#TODO MAKE "FORK PLOTS of mutation prevalence by number of mutations


print len(endometrialHypermutators)











