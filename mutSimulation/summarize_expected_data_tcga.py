#written by Noah Friedman
#basically the same as summarize_expected_data.py except it is used for TCGA data and is thus formatted differently
import sys
import argparse
import os
import pandas as pd
import numpy as np
import re

sys.path.append('/ifs/work/taylorlab/friedman/myUtils')
import analysis_utils 
import mutationSigUtils
import mutation_modeling_util

def calculate_oncogenic_mut_susceptibility_of_genes_by_signature(oncogenicSDict, suffix='_hotspot_rate'):
    listOfDicts = []
    sigNames = ['Signature.' + str(i) for i in range(1,31)]
    for i in range(1,31):
        curSig = 'Signature.' + str(i)
        d = {}
        for s in sigNames:
            d[s] = 0
        d[curSig] = 1
        #PRETEND we got a case with 100% signature i on the decomposition
        quadNucFractions = mutation_modeling_util.get_quadnuc_fracs_given_decomposition(d, spectraPath = pathPrefix + '/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt')
        #v = mutation_modeling_util.get_expected_oncogenic_val_given_quadnuc_fractions(quadNucFractions, oncogenicSDict, 'IMPACT_468')
        #ALERT NOAH I CHANGED THIS HERE
        v = mutation_modeling_util.get_expected_oncogenic_val_given_quadnuc_fractions_v2(quadNucFractions, oncogenicSDict, suffix)

        listOfDicts.append({'Signature_Name': curSig, 'ExpectedFracOfMutsOncogenic': v})
    return pd.DataFrame(listOfDicts)

def quantify_quadnuc_oncogenic_susceptibility_per_mutation(simMafData, relatedGenes=None):
    allBases = ['A', 'C', 'G', 'T']
    changes = ['CA', 'CG', 'CT', 'TA', 'TC', 'TG'] #format: 'CA' means a change from C>A
    allQuadNucs = [firstBase + change + lastBase for firstBase in allBases for change in changes for lastBase in allBases] #enumerate all 96 quadnucs for signatures
   
    d = {}
    
    for quadNuc in allQuadNucs:
        nPossibleMuts = sum(simMafData[quadNuc]) - sum(simMafData[quadNuc + '_silent'])
        nPossibleHotspotMuts = sum(simMafData[quadNuc + '_oncogenic'])
        d[quadNuc + '_oncogenic_rate'] = (1.0*nPossibleHotspotMuts)/(1.0*nPossibleMuts)
        
    return d



tcgaSigsPath = '/ifs/work/taylorlab/pensona/dmp_sigs/tcga/mc3.v0.2.8.PUBLIC.LAML_PATCH.trinuc.30sigs.txt'
tcgaSigs = pd.read_table(tcgaSigsPath)
tcgaSigs['sampleNameAdj'] = tcgaSigs['Sample Name'].apply(lambda x: x[:15])

#load in nmut info
nmutDict = {}
nmutLines = open('/ifs/work/taylorlab/friedman/myAdjustedDataFiles/tcgaNonSynomCounts.tsv').readlines()
for line in nmutLines:
	
	if len(line.strip('\n').split('\t')) == 2:
		case, nmut = line.strip('\n').split('\t')
		nmutDict[case] = float(nmut)

tcgaSigs['nNonSynom'] = tcgaSigs['sampleNameAdj'].apply(lambda x: nmutDict[x] if x in nmutDict else None)

tcgaSigs['Signature.MMR'] = tcgaSigs['Signature.6'] + tcgaSigs['Signature.15'] + tcgaSigs['Signature.20']+ tcgaSigs['Signature.21'] + tcgaSigs['Signature.26']
tcgaSigs['Signature.POLE'] = tcgaSigs['Signature.10'] + tcgaSigs['Signature.14']
tcgaSigs['mmrAttributed'] = tcgaSigs['Signature.MMR'] * tcgaSigs['nNonSynom']
tcgaSigs['poleAttributed'] = tcgaSigs['Signature.POLE'] * tcgaSigs['nNonSynom']
mmrCases = set(tcgaSigs[tcgaSigs['mmrAttributed'] > 150]['sampleNameAdj'])
poleCases = set(tcgaSigs[tcgaSigs['poleAttributed'] > 150]['sampleNameAdj'])

hyperCases = mmrCases | mmrCases

print 'loading spectra d'
spectraPath = '/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt'
spectraD = mutationSigUtils.convert_spectrum_file_to_dict_of_dicts(spectrumFile=spectraPath)

print 'loading expected mut rate df'
expectedMutRateDf = pd.read_table('/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/mutSimulation/expectedMutationTables/quadNucAndGenePossibleMutationSummaryIMPACT341.tsv')
genes = set(expectedMutRateDf['gene'])
#save runtime by presubsetting the dataframe by genes
gdfDict = {}
for gene in genes:
	gdf = expectedMutRateDf[expectedMutRateDf['gene'] == gene]
	gdfDict[gene] = gdf

listOfDfs = []
cntr = 0

print 'beginning iteration over cases, note each iteration takes between 5-10 and seconds per case'
for case in list(hyperCases):
	print 'analyzing case number ', cntr, ' out of ', len(hyperCases), ' cases'
	cntr += 1
	caseSigs = tcgaSigs[tcgaSigs['sampleNameAdj'] == case]
	quadNucFracs = mutation_modeling_util.get_quadnuc_fracs_given_decomposition(caseSigs.iloc[0], spectraD)

	expectedCopy = expectedMutRateDf.copy()
	expectedCopy['hotspotChance'] =  expectedCopy.apply(lambda row: quadNucFracs[row['quadNuc']]* row['hotspotChance'], axis=1)
	expectedCopy['oncogenicChance'] =  expectedCopy.apply(lambda row: quadNucFracs[row['quadNuc']]* row['oncogenicChance'], axis=1)
	expectedCopy['truncatingChance'] =  expectedCopy.apply(lambda row: quadNucFracs[row['quadNuc']]* row['truncatingChance'], axis=1)

	#THE Groupby operation gives us the information we want
	hotspotInfo = expectedCopy.groupby("gene").hotspotChance.sum().reset_index()
	oncogenicInfo = expectedCopy.groupby("gene").oncogenicChance.sum().reset_index()
	truncatingInfo = expectedCopy.groupby("gene").truncatingChance.sum().reset_index()

	#merge them all
	dfMerged = reduce(lambda left,right: pd.merge(left,right,on='gene'), [hotspotInfo, oncogenicInfo, truncatingInfo])
	dfMerged['case'] = case
	listOfDfs.append(dfMerged)

df = pd.concat(listOfDfs)
print 'writing data'
df.to_csv('/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/mutSimulation/expectedMutationTables/tcgaHypermutatorsExpectedGeneMutInfo.tsv', index=False, sep='\t')













