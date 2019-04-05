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

def calculate_genomic_region_impact_size_info(chrInfo, cdsInfo, cdsIntervals):
    
    chrInfo['arm'] = chrInfo['Location'].apply(lambda x: x.split('p')[0] + 'p' if 'p' in x else x.split('q')[0] + 'q')
    armDict = dict(zip(chrInfo['Gene_Symbol'], chrInfo['arm']))
    cdsDict = dict(zip(cdsInfo['Hugo_Symbol'], cdsInfo['cds_length']))
    
    armPanelContentDict = {}
    for gene in set(cdsInfo['Hugo_Symbol']):
       if gene in armDict:
            arm = armDict[gene]
            if arm not in armPanelContentDict:
                armPanelContentDict[arm] = cdsDict[gene]
            else:
                armPanelContentDict[arm] = armPanelContentDict[arm] + cdsDict[gene]
    return armPanelContentDict, armDict

def asses_case_per_region_mut_info(muts, armContentDict):
    listOfDicts = []
    for case in set(muts['Tumor_Sample_Barcode']):
        caseMuts = muts[muts['Tumor_Sample_Barcode'] == case]
        for region in set(caseMuts['chrArm']):
            if region in armContentDict:
                listOfDicts.append({
                        'mutRate': 1.0*caseMuts[caseMuts['chrArm'] == region].shape[0]/armContentDict[region],
                        'mutsInRegion': 1.0*caseMuts[caseMuts['chrArm'] == region].shape[0],
                        'region': region,
                        'case': case
                })
    return pd.DataFrame(listOfDicts)


chromosomePosInfo = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/IMPACT_genes.txt')
#cytobandInfo = cytobandInfo.rename(columns={'X0': 'Chromosome', 'X1': 'Start_Pos', 'X2': 'End_Pos', 'X3': 'ChromPosition', 'X4': 'X4'})


cdsInfo = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/impact_gene_reference_signatures.tsv')
cdsIntervals = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/genelist.with_aa.interval_list', header=None, prefix='X')
cdsIntervals = cdsIntervals.rename(columns={'X0': 'Chromosome', 'X1': 'Start_Pos', 'X2': 'End_Pos', 'X3': 'Strand', 'X4': 'LabelAndExon'})
cdsIntervals['Hugo_Symbol'] = cdsIntervals['LabelAndExon'].apply(lambda x: str(x).split(':')[0])

armContentDict, armDictInfo = calculate_genomic_region_impact_size_info(chromosomePosInfo, cdsInfo, cdsIntervals)
armDictInfo['KMT2A'] = armDictInfo['MLL']
armDictInfo['KMT2D'] = armDictInfo['MLL2']
armDictInfo['KMT2C'] = armDictInfo['MLL3']

####AFTER initiating all the arm content info please go to looking at real data

mutationMaf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/hypermutationAnalysisProjectData/hypermutationAnalysisProjectMaf.tsv')
impactSigs = pd.read_table(pathPrefix + '/ifs/res/taylorlab/impact_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
impactSigs['pid'] = impactSigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
impactSigs['cancer_type'] = impactSigs['pid'].apply(lambda x: cDict[x] if x in cDict else None)

cntr = 0
for x in set(impactSigs[impactSigs['cancer_type'] == 'Colorectal Cancer']['Tumor_Sample_Barcode']):
    cntr += 1
    if cntr%50 == 0: print x
    
mutationMaf['oncogenic']
    
for i in set(mutationMaf[(mutationMaf['Hugo_Symbol'] == 'POLE') & (mutationMaf['oncogenic'].notnull())]['Tumor_Sample_Barcode']):
    print i
    
print impactSigs[impactSigs['seqType'] == 'IH3']
impactSigs['seqType'] = impactSigs['Tumor_Sample_Barcode'].apply(lambda x: x[-3:])    


gliomaHypermutants = set(impactSigs[(impactSigs['cancer_type'] == 'Glioma') & (impactSigs['Nmut_Mb'] > 75)]['Tumor_Sample_Barcode'])
gliomaHypermutantMutations = mutationMaf[mutationMaf['Tumor_Sample_Barcode'].isin(gliomaHypermutants)]

gliomaHypermutantMutations['chrArm'] = gliomaHypermutantMutations['Hugo_Symbol'].apply(lambda x: armDictInfo[x] if x in armDictInfo else None)


df = asses_case_per_region_mut_info(gliomaHypermutantMutations, armContentDict)
df['orderingVal'] = df['region'].apply(lambda x: 
    23 if x == 'Xp' else 23.1 if x == 'Xq'
    else int(x.split('p')[0]) if 'p' in x else int(x.split('q')[0]))
df.to_csv('~/Desktop/WORK/dataForLocalPlotting/locationDist.tsv', sep='\t', index=False)




