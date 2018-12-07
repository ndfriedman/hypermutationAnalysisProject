#written by Noah Friedman (a template for scripts to be excuted in the spyder environment
#ALERT THIS IS THE OFFLINE VERSION FYI
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

#returns a DF counting the number of TOP5, 5-10, and 10-20 tmz genes
def enumerate_pre_and_post_glioma_oncogenic_muts(preTMZMaf, postTMZMaf, top5TMZ, fiveThrough10TMZ, tenThrough20TMZ):
    listOfDicts = []
    top20TMZGenes = top5TMZ | fiveThrough10TMZ | tenThrough20TMZ
    
    oncoKbOncogenicAnnotations = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
    preTMZoncogenicMuts = preTMZMaf[preTMZMaf['oncogenic'].isin(oncoKbOncogenicAnnotations)]
    postTMZOncogenicMuts = postTMZMaf[postTMZMaf['oncogenic'].isin(oncoKbOncogenicAnnotations)]
    preTMZgenes = set(preTMZoncogenicMuts['Hugo_Symbol'])
    postTMZgenes = set(postTMZOncogenicMuts['Hugo_Symbol']) 
    bothSamplesGenes = preTMZgenes & postTMZgenes
    onlyPostTMZGenes = postTMZgenes - preTMZgenes
    
    #SET up long format DF entries 
    tsb = preTMZMaf['Tumor_Sample_Barcode'].iloc[0]
    listOfDicts.append({'Tumor_Sample_Barcode': tsb, 'NPreTMZOncogenicMuts': len(bothSamplesGenes & top5TMZ), 'NPostTMZOncogenicMuts': len(onlyPostTMZGenes & top5TMZ), 'class': 'Top 5 Most Mutated Genes In Glioma'})
    listOfDicts.append({'Tumor_Sample_Barcode': tsb, 'NPreTMZOncogenicMuts': len(bothSamplesGenes & fiveThrough10TMZ), 'NPostTMZOncogenicMuts': len(onlyPostTMZGenes & fiveThrough10TMZ), 'class': '6th-10th Most Mutated Genes In Glioma'})
    listOfDicts.append({'Tumor_Sample_Barcode': tsb, 'NPreTMZOncogenicMuts': len(bothSamplesGenes & tenThrough20TMZ), 'NPostTMZOncogenicMuts': len(onlyPostTMZGenes & tenThrough20TMZ), 'class': '11th-20th Most Mutated Genes In Glioma'})
    listOfDicts.append({'Tumor_Sample_Barcode': tsb, 'NPreTMZOncogenicMuts': len(bothSamplesGenes - top20TMZGenes), 'NPostTMZOncogenicMuts': len(onlyPostTMZGenes - top20TMZGenes), 'class': 'Genes infrequently mutated in Glioma'})
    
    return listOfDicts #set up for a long format DF

#TODO fix issues with cases that have 2X post or pre samples
#starts out an embryonic df with tumor sample barcode, gene, and a pre/post tmz mutaation flag
    
#TODO--also show genes that are only detected in the pre tmz MAF
def enumerate_pre_and_post_glioma_oncogenic_muts_mode_2(preTMZMaf, postTMZMaf):
    listOfDicts = []
    
    oncoKbOncogenicAnnotations = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
    preTMZoncogenicMuts = preTMZMaf[preTMZMaf['oncogenic'].isin(oncoKbOncogenicAnnotations)]
    postTMZOncogenicMuts = postTMZMaf[postTMZMaf['oncogenic'].isin(oncoKbOncogenicAnnotations)]
    preTMZgenes = set(preTMZoncogenicMuts['Hugo_Symbol'])
    postTMZgenes = set(postTMZOncogenicMuts['Hugo_Symbol']) 
    bothSamplesGenes = preTMZgenes & postTMZgenes
    onlyPostTMZGenes = postTMZgenes - preTMZgenes
    
    #SET up long format DF entries 
    tsb = preTMZMaf['Tumor_Sample_Barcode'].iloc[0]
    for gene in bothSamplesGenes:
        listOfDicts.append({'Tumor_Sample_Barcode': tsb,  'Gene': gene, 'preTMZ': True, 'postTMZ': False})
    for gene in onlyPostTMZGenes:
        listOfDicts.append({'Tumor_Sample_Barcode': tsb,  'Gene': gene, 'preTMZ': False, 'postTMZ': True})
    return listOfDicts #set up for a long format DF

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
        caseDfUncommonPre = caseDf[caseDf['uncommonGliomaGeneDisplayNamePreTMZ'].notnull()]
        caseDfUncommonPost = caseDf[caseDf['uncommonGliomaGeneDisplayNamePostTMZ'].notnull()]
        
        preTMZOrder = enumerate_order(set(caseDfUncommonPre['Gene']), rankingDict)
        for entry in preTMZOrder.keys():
            dictKey = case + '_' + entry
            lambdaDict[dictKey] = preTMZOrder[entry]

        postTMZOrder = enumerate_order(set(caseDfUncommonPost['Gene']), rankingDict)
        for entry in postTMZOrder.keys():
            dictKey = case + '_' + entry
            lambdaDict[dictKey] = postTMZOrder[entry]
    
    return lambdaDict

#defines a way of ordering cases for an oncoprint type plot by looking at the topNGenes and seeing if they are in there
#cases are ordered by preTMZ mut status
def define_case_ordering(df, topNGenes):
    cases = set(df['Tumor_Sample_Barcode'])
    df = df[df['preTMZ'] == True] #only take preTMZ mutated genes
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
            
        

#a little function to help me out finding out stuff
def create_lambda_dict_key(row, colName):
    if row['Tumor_Sample_Barcode'] == None or row[colName] == None:
        return None
    else:
        return row['Tumor_Sample_Barcode'] + '_' + row[colName]
        
#a function whose existence is very annoying because I struggle with ggplot
#disgusting code fyi; make better plz
        
#(xAxisDisplay, -fracMutated), y=reorder(Tumor_Sample_Barcode*, caseRank*), fill=fracMutated*  textDisplay*    
    
def fix_uncommon_col_for_ggplot(df):
    dfsToAdd = []
    cases = set(df['Tumor_Sample_Barcode'])
    for case in cases:
        caseDf = df[df['Tumor_Sample_Barcode'] == case]
        if caseDf[caseDf['uncommonGliomaGeneDisplayNamePreTMZ'].notnull()].shape[0] == 0:
            #we add a single dummy row to the dataframe 
            caseRank = caseDf.iloc[0]['caseRank']
            dummyDf = pd.DataFrame([{'Tumor_Sample_Barcode': case, 'caseRank': caseRank, 'fracMutated':0, 'preTMZUncommonOrder': 0, 'Gene':''}])
            dfsToAdd.append(dummyDf)
    return pd.concat([df]+ dfsToAdd)

mafWithInfo = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/annotatedOncoPlusHotspotMafAllImpact_trinuc')
mafWithInfo['pid'] = mafWithInfo['Tumor_Sample_Barcode'].apply(lambda x: x[:9])

impactSigs = pd.read_table(pathPrefix + '/ifs/res/taylorlab/impact_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
impactSigs['pid'] = impactSigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
impactSigs['cancer_type'] = impactSigs['pid'].apply(lambda x: cDict[x] if x in cDict else None)

gliomaSigs = impactSigs[impactSigs['cancer_type'] == 'Glioma']
gliomasWithTMZ = gliomaSigs[(gliomaSigs['mean_11'] >= .25) & (gliomaSigs['Nmut'] >= 10)]
gliomasWithoutTMZ = gliomaSigs[(gliomaSigs['mean_11'] < .25) | (gliomaSigs['Nmut'] <= 10)]
tmzEvidentPatients = set(gliomasWithTMZ['pid'])

twoSamplePreTMZGliomas = gliomasWithoutTMZ[gliomasWithoutTMZ['pid'].isin(tmzEvidentPatients)]
twoSamplePatients = set(twoSamplePreTMZGliomas['pid'])
twoSamplePostTMZGliomas = gliomasWithTMZ[gliomasWithTMZ['pid'].isin(twoSamplePatients)]

gliomaMuts = mafWithInfo[mafWithInfo['Tumor_Sample_Barcode'].isin(set(gliomaSigs['Tumor_Sample_Barcode']))]

print gliomaMuts[gliomaMuts['Tumor_Sample_Barcode'] == 'P-0000883-T01-IM3']['oncogenic']


nGenesMutatedInGlioma = len(set(gliomaMuts['Hugo_Symbol']))

#divide every gene in the fractional dict by the n glioma cases
gliomaGenesFractionalDict = {key: value for (key, value) in topGliomaGenes}
for key, value in gliomaGenesFractionalDict.items():
    gliomaGenesFractionalDict[key] = 1.0*value/gliomaSigs.shape[0]

topGliomaGenes = maf_analysis_utils.enumerate_top_n_oncogenic_mutated_genes_across_cohort(gliomaMuts, nGenesMutatedInGlioma)
topGliomaGeneNames = set(['TERT', 'TP53', 'IDH1', 'PTEN', 'ATRX'])

#make a dictionary of rankings
gliomaGeneRanking = dict()
cntr = 1
for gene in topGliomaGenes:
    gliomaGeneRanking[gene[0]] = cntr
    cntr += 1

listOfDicts = []
for patient in twoSamplePatients:
    preIds = set(twoSamplePreTMZGliomas[twoSamplePreTMZGliomas['pid'] == patient]['Tumor_Sample_Barcode'])
    postIds = set(twoSamplePostTMZGliomas[twoSamplePostTMZGliomas['pid'] == patient]['Tumor_Sample_Barcode'])
    for preId in preIds:
        for postId in postIds:
            preMuts = gliomaMuts[gliomaMuts['Tumor_Sample_Barcode'] == preId]
            postMuts = gliomaMuts[gliomaMuts['Tumor_Sample_Barcode'] == postId]
            #listOfDicts += enumerate_pre_and_post_glioma_oncogenic_muts(preMuts, postMuts, top5TMZGenes, fiveThrough10TMZGenes, tenThrough20TMZGenes)
            listOfDicts += enumerate_pre_and_post_glioma_oncogenic_muts_mode_2(preMuts, postMuts)

df = pd.DataFrame(listOfDicts)
df['postTMZGene'] = df.apply(lambda row: row['Gene'] if row['postTMZ'] == True else None, axis=1)
df['preTMZGene'] = df.apply(lambda row: row['Gene'] if row['preTMZ'] == True else None, axis=1)
df['topGliomaGeneDisplayNamePreTMZ'] = df.apply(lambda row: row['Gene'] if (row['Gene'] in topGliomaGeneNames and row['preTMZ'] == True) else None, axis = 1)
df['topGliomaGeneDisplayNamePostTMZ'] = df.apply(lambda row: row['Gene'] if (row['Gene'] in topGliomaGeneNames and row['postTMZ'] == True) else None, axis = 1)

df['uncommonGliomaGeneDisplayNamePreTMZ'] = df.apply(lambda row: row['Gene'] if (row['Gene'] not in topGliomaGeneNames and row['preTMZ'] == True) else None, axis = 1)
df['uncommonGliomaGeneDisplayNamePostTMZ'] = df.apply(lambda row: row['Gene'] if (row['Gene'] not in topGliomaGeneNames and row['postTMZ'] == True) else None, axis = 1)

df['fracMutated'] = df['Gene'].apply(lambda x: gliomaGenesFractionalDict[x] if x in gliomaGenesFractionalDict else 0)

lDict = add_case_gene_mut_ranking(df, gliomaGeneRanking)

df['preTMZUncommonOrder'] = df.apply(lambda row: lDict[create_lambda_dict_key(row, 'uncommonGliomaGeneDisplayNamePreTMZ')] if create_lambda_dict_key(row, 'uncommonGliomaGeneDisplayNamePreTMZ') in lDict else None, axis=1)
df['postTMZUncommonOrder'] = df.apply(lambda row: lDict[create_lambda_dict_key(row, 'uncommonGliomaGeneDisplayNamePostTMZ')] if create_lambda_dict_key(row, 'uncommonGliomaGeneDisplayNamePostTMZ') in lDict else None, axis=1)

top5Ranking = define_case_ordering(df, ['TERT', 'TP53', 'IDH1', 'PTEN', 'ATRX'])
df['caseRank'] = df['Tumor_Sample_Barcode'].apply(lambda x: top5Ranking[x])
#row['Tumor_Sample_Barcode'] + '_' + row['uncommonGliomaGeneDisplayNamePreTMZ']

df = fix_uncommon_col_for_ggplot(df)


df.to_csv('~/Desktop/dataForLocalPlotting/tmzMultipleSampleGenes.tsv', sep='\t', index=False)

#TODO GLIOMA pre vs POST Analysis of cases before and after recurrence



