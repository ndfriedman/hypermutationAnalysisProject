#written by Noah Friedman (a template for scripts to be excuted in the spyder environment
#try to compare distributions of oncogenic mutations in hypermutators vs not
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

#given motifGroups which is a 
def segement_mutations_by_motif(motifGroups):
    return 0

def do_glioma_TMZ_analysis(mafWithInfo, impactSigs):
    gliomaSigs = impactSigs[impactSigs['cancer_type'] == 'Glioma']
    gliomasWithTMZ = gliomaSigs[(gliomaSigs['mean_11'] >= .25) & (gliomaSigs['Nmut'] >= 10)]
    gliomasWithTMZMuts = mafWithInfo[mafWithInfo['Tumor_Sample_Barcode'].isin(set(gliomasWithTMZ['Tumor_Sample_Barcode']))]
    gliomasWithTMZMuts['quadNuc'] = gliomasWithTMZMuts.apply(lambda row: mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)
    
    allGliomaMuts = mafWithInfo[mafWithInfo['Tumor_Sample_Barcode'].isin(set(gliomaSigs['Tumor_Sample_Barcode']))]

    tmzTopMotifs = set(['ACTC', 'CCTC', 'GCTC', 'TCTC'])
    tmzSecondTierMotifs = set(['ACTT', 'CCTT', 'GCTT', 'TCTT'])
    tmzThirdTierMotifs = set(['ACTA', 'CCTA', 'GCTA', 'TCTA'])
    allOtherMotifs = analysis_utils.get_all_possible_quadNucs() - tmzTopMotifs - tmzSecondTierMotifs - tmzThirdTierMotifs
    
    oncoKbOncogenicAnnotations = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
    tmzOncogenicMuts = gliomasWithTMZMuts[gliomasWithTMZMuts['oncogenic'].isin(oncoKbOncogenicAnnotations)]
    
    mutationsAttmzTopMotifs = tmzOncogenicMuts[tmzOncogenicMuts['quadNuc'].isin(tmzTopMotifs)]
    mutationsAttmzSecondTierMotifs = tmzOncogenicMuts[tmzOncogenicMuts['quadNuc'].isin(tmzSecondTierMotifs)]
    mutationsAttmzThirdTierMotifs = tmzOncogenicMuts[tmzOncogenicMuts['quadNuc'].isin(tmzThirdTierMotifs)]
    mutationsAtAllOtherMotifs = tmzOncogenicMuts[tmzOncogenicMuts['quadNuc'].isin(allOtherMotifs)]
    
    listOfDicts = []
    for index, row in mutationsAttmzTopMotifs.iterrows():
        listOfDicts.append({'Tumor_Sample_Barcode': row['Tumor_Sample_Barcode'], 'Gene': row['Hugo_Symbol'], 'motif': 'TMZTop'})
    for index, row in mutationsAttmzSecondTierMotifs.iterrows():
        listOfDicts.append({'Tumor_Sample_Barcode': row['Tumor_Sample_Barcode'], 'Gene': row['Hugo_Symbol'], 'motif': 'TMZSecond'})
    for index, row in mutationsAttmzThirdTierMotifs.iterrows():
        listOfDicts.append({'Tumor_Sample_Barcode': row['Tumor_Sample_Barcode'], 'Gene': row['Hugo_Symbol'], 'motif': 'TMZThird'})
    for index, row in mutationsAtAllOtherMotifs.iterrows():
        listOfDicts.append({'Tumor_Sample_Barcode': row['Tumor_Sample_Barcode'], 'Gene': row['Hugo_Symbol'], 'motif': 'Other'})
    
    return pd.DataFrame(listOfDicts)


def do_pole_analysis(mafWithInfo, poleSigs, poleCaseMuts):
    
    spectraEnrichmentDict = mutationSigUtils.get_enriched_spectra_for_signatures(spectraSignificanceThresh=.05, pathPrefix='/Users/friedman/Desktop/mnt',
	signaturesToIgnore= #ignore signatures we dont care about 
	set(['Signature.5','Signature.8','Signature.9','Signature.12','Signature.16','Signature.19','Signature.22','Signature.23','Signature.24','Signature.25','Signature.27','Signature.28','Signature.29','Signature.30']))

    
    poleSpectra = spectraEnrichmentDict['Signature.10']
    sig14Spectra = spectraEnrichmentDict['Signature.14']
    agingSpectra = spectraEnrichmentDict['Signature.1']
    
    #these are the four categories we will annotate stuff for
    poleAndAgingSpectra = poleSpectra & agingSpectra
    agingOnlySpectra = agingSpectra - poleSpectra
    poleOnlySpectra = poleSpectra - agingSpectra
    sig14OnlySpectra = sig14Spectra - poleSpectra - agingSpectra
    otherSpectra = analysis_utils.get_all_possible_quadNucs() - poleAndAgingSpectra - agingOnlySpectra - poleOnlySpectra - sig14OnlySpectra
    
    oncoKbOncogenicAnnotations = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
    poleOncogenicMuts = poleCaseMuts[poleCaseMuts['oncogenic'].isin(oncoKbOncogenicAnnotations)]
     
    poleAndAgingMotifMuts = poleOncogenicMuts[poleOncogenicMuts['quadNuc'].isin(poleAndAgingSpectra)]
    agingOnlyMotifMuts = poleOncogenicMuts[poleOncogenicMuts['quadNuc'].isin(agingOnlySpectra)]
    poleOnlyMotifMuts = poleOncogenicMuts[poleOncogenicMuts['quadNuc'].isin(poleOnlySpectra)]
    sig14OnlyMotifMuts = poleOncogenicMuts[poleOncogenicMuts['quadNuc'].isin(sig14OnlySpectra)]
    otherMotifMuts = poleOncogenicMuts[poleOncogenicMuts['quadNuc'].isin(otherSpectra)]
    
    listOfDicts = []
    for index, row in poleAndAgingMotifMuts.iterrows():
        listOfDicts.append({'Tumor_Sample_Barcode': row['Tumor_Sample_Barcode'], 'Gene': row['Hugo_Symbol'], 'motif': 'PoleAndAging', 'Hotspot': row['is-a-hotspot']})
    for index, row in agingOnlyMotifMuts.iterrows():
        listOfDicts.append({'Tumor_Sample_Barcode': row['Tumor_Sample_Barcode'], 'Gene': row['Hugo_Symbol'], 'motif': 'AgingOnly', 'Hotspot': row['is-a-hotspot']})
    for index, row in poleOnlyMotifMuts.iterrows():
        listOfDicts.append({'Tumor_Sample_Barcode': row['Tumor_Sample_Barcode'], 'Gene': row['Hugo_Symbol'], 'motif': 'PoleOnly', 'Hotspot': row['is-a-hotspot']})
    for index, row in sig14OnlyMotifMuts.iterrows():
        listOfDicts.append({'Tumor_Sample_Barcode': row['Tumor_Sample_Barcode'], 'Gene': row['Hugo_Symbol'], 'motif': 'Sig14Only', 'Hotspot': row['is-a-hotspot']})
    for index, row in otherMotifMuts.iterrows():
        listOfDicts.append({'Tumor_Sample_Barcode': row['Tumor_Sample_Barcode'], 'Gene': row['Hugo_Symbol'], 'motif': 'Other', 'Hotspot': row['is-a-hotspot']})
    
    return pd.DataFrame(listOfDicts)
    
    
    
def assign_plot_rank(row, rankDict):
    if row['Gene'] in rankDict:
        return rankDict[row['Gene']]
    else:
        return len(rankDict) + 1
        
#defines a way of ordering cases for an oncoprint type plot by looking at the topNGenes and seeing if they are in there
#cases are ordered by preTMZ mut status
def define_case_ordering(df, topNGenes):
    cases = set(df['Tumor_Sample_Barcode'])
    rankingDict = dict()
    for case in cases:
        rankingDict[case] = 0
    for case in cases:
        localMuts = df[df['Tumor_Sample_Barcode'] == case]
        cntr = 0
        caseRank = 0
        for gene in topNGenes:
            if localMuts[localMuts['Gene'] == gene].shape[0] > 0:
                caseRank += 10.0**(-1*cntr)
            cntr += 1
        rankingDict[case] = caseRank
    return rankingDict

def add_case_gene_mut_ranking(df, rankingDict):
    
    lambdaDict = dict() #HERES a dict that we return for later use as a lambda function to make stuff work
    def enumerate_order(vals, rDict):
        l = [(v, rDict[v]) for v in vals]
        sortedL = sorted(l, key=lambda x: x[1])
        returnD = dict()
        cntr = 0
        for v in sortedL:
            returnD[v[0]] =  cntr
            cntr += 1
        return returnD
        
    cases = set(df['Tumor_Sample_Barcode'])
    for case in cases:
        caseDf = df[df['Tumor_Sample_Barcode'] == case]
        
        order = enumerate_order(set(caseDf['Gene']), rankingDict)
        for entry in order.keys():
            dictKey = case + '_' + entry
            lambdaDict[dictKey] = order[entry]

    return lambdaDict


#a little function to help me out finding out stuff
def create_lambda_dict_key(row, colName):
    if row['Tumor_Sample_Barcode'] == None or row[colName] == None:
        return None
    else:
        return row['Tumor_Sample_Barcode'] + '_' + row[colName]

            
        

mafWithInfo = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/annotatedOncoPlusHotspotMafAllImpact_trinuc')
mafWithInfo['pid'] = mafWithInfo['Tumor_Sample_Barcode'].apply(lambda x: x[:9])

impactSigs = pd.read_table(pathPrefix + '/ifs/res/taylorlab/impact_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
impactSigs['pid'] = impactSigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
impactSigs['cancer_type'] = impactSigs['pid'].apply(lambda x: cDict[x] if x in cDict else None)


#GLIOMA Analysis
occurenceCounter, rankingDict, gliomaGenesFractionalDict = maf_analysis_utils.enumerate_top_n_oncogenic_mutated_genes_across_cohort(allGliomaMuts)
gliomaMutTimingDf = do_glioma_TMZ_analysis(mafWithInfo, impactSigs)
gliomaMutTimingDf['ranking'] = gliomaMutTimingDf['Gene'].apply(lambda x: rankingDict[x])
top5GliomaGenes = ['TERT', 'TP53', 'IDH1', 'PTEN', 'ATRX']
gliomaMutTimingDf['fracMutated'] = gliomaMutTimingDf['Gene'].apply(lambda x: gliomaGenesFractionalDict[x] if x in gliomaGenesFractionalDict else 0)
top5Ranking = define_case_ordering(gliomaMutTimingDf, ['TERT', 'TP53', 'IDH1', 'PTEN', 'ATRX'])
gliomaMutTimingDf['caseRank'] = gliomaMutTimingDf['Tumor_Sample_Barcode'].apply(lambda x: top5Ranking[x])
gliomaMutTimingDf['topGeneDisplayName'] = gliomaMutTimingDf['Gene'].apply(lambda x: x if x in top5GliomaGenes else None)
gliomaMutTimingDf['infrequentGeneDisplayName'] = gliomaMutTimingDf['Gene'].apply(lambda x: x if x  not in top5GliomaGenes else None)
lDict = add_case_gene_mut_ranking(gliomaMutTimingDf[gliomaMutTimingDf['infrequentGeneDisplayName'].notnull()], rankingDict)
gliomaMutTimingDf['uncommonOrder'] = gliomaMutTimingDf.apply(lambda row: lDict[create_lambda_dict_key(row, 'infrequentGeneDisplayName')] if create_lambda_dict_key(row, 'infrequentGeneDisplayName') in lDict else None, axis=1)
gliomaMutTimingDf.to_csv('~/Desktop/dataForLocalPlotting/tmzCasesByMotif.tsv', sep='\t', index=False)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ENDOMETRIAL POLE ANALYSIS
endometrialSigs = impactSigs[impactSigs['cancer_type'] == 'Endometrial Cancer']
endometrialPoleSigs = endometrialSigs[(endometrialSigs['mean_10'] > .15) & (endometrialSigs['Nmut'] > 100)]
allEndometrialMuts = mafWithInfo[mafWithInfo['Tumor_Sample_Barcode'].isin(set(endometrialSigs['Tumor_Sample_Barcode']))]
endometrialPoleCaseMuts = allEndometrialMuts[allEndometrialMuts['Tumor_Sample_Barcode'].isin(set(endometrialPoleSigs['Tumor_Sample_Barcode']))]
endometrialPoleCaseMuts['quadNuc'] = endometrialPoleCaseMuts.apply(lambda row: mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)
occurenceCounter, rankingDict, endometrialGenesFractionalDict = maf_analysis_utils.enumerate_top_n_oncogenic_mutated_genes_across_cohort(allEndometrialMuts)

endometrialPoleMutDataDf = do_pole_analysis(allEndometrialMuts, endometrialPoleSigs, endometrialPoleCaseMuts)
top10endometrialMuts = [x[0] for x in Counter(endometrialGenesFractionalDict).most_common(10)]
top10endometrialMuts = ['POLE'] + top10endometrialMuts
rankingDict['POLE'] = 0

endometrialPoleMutDataDf['ranking'] = endometrialPoleMutDataDf['Gene'].apply(lambda x: rankingDict[x])
endometrialPoleMutDataDf['fracMutated'] = endometrialPoleMutDataDf['Gene'].apply(lambda x: endometrialGenesFractionalDict[x] if x in endometrialGenesFractionalDict else 0)
top10Ranking = define_case_ordering(endometrialPoleMutDataDf, top10endometrialMuts)
endometrialPoleMutDataDf['caseRank'] = endometrialPoleMutDataDf['Tumor_Sample_Barcode'].apply(lambda x: top10Ranking[x])
endometrialPoleMutDataDf['topGeneDisplayName'] = endometrialPoleMutDataDf['Gene'].apply(lambda x: x if x in top10endometrialMuts else None)
endometrialPoleMutDataDf['infrequentGeneDisplayName'] = endometrialPoleMutDataDf['Gene'].apply(lambda x: x if x  not in top10endometrialMuts else None)

lDict = add_case_gene_mut_ranking(endometrialPoleMutDataDf[endometrialPoleMutDataDf['infrequentGeneDisplayName'].notnull()], rankingDict)
endometrialPoleMutDataDf['uncommonOrder'] = endometrialPoleMutDataDf.apply(lambda row: lDict[create_lambda_dict_key(row, 'infrequentGeneDisplayName')] if create_lambda_dict_key(row, 'infrequentGeneDisplayName') in lDict else None, axis=1)
endometrialPoleMutDataDf.to_csv('~/Desktop/dataForLocalPlotting/poleEndometrialCasesByMotif.tsv', sep='\t', index=False)



Counter(endometrialGenesFractionalDict).most_common(10)

