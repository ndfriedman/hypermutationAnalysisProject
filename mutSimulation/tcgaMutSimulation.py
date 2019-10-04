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
import mutation_modeling_util

tcgaMc3Maf = pd.read_table(pathPrefix + '/ifs/res/taylorlab/ang46/ext/mafs/mc3/mc3.v0.2.8.PUBLIC.LAML_PATCH_prepped_facets_oncokb.maf')

tcgaMc3Maf.columns.values

caseDecompositionDict = {}

signatureCols = ['Signature.' + str(i) for i in range(1,31)]
for col in signatureCols:
    caseDecompositionDict[col] = 0
    
caseDecompositionDict['Signature.10'] = 1

ucecSignatures = pd.read_table(pathPrefix + '/ifs/res/taylorlab/jonssonp/msk_impact_gml_som/tcga/mutsig/UCEC.mc3.v0.2.8.PUBLIC.maf.2018-01-29.oncokb.vep.mutsig.txt')
ucecSignatures['betterId'] = ucecSignatures['Sample Name'].apply(lambda x: x[:15])
cases = ucecSignatures['Sample Name']
ucecSignatures.columns.values

mutFracDictVals = mutation_modeling_util.get_simulated_oncogenic_mut_frac_per_case(
        ucecSignatures, 
        spectraD, 
        quadNucDict,
        quadNucDictAdjusted,
        geneOncogenicOrHotspotProbDict,
        geneMutProbInfo,
        n=10000,                                                         
        idColName='Sample Name')

ucecMuts = tcgaMc3Maf[tcgaMc3Maf['SAMPLE_ID'].isin(set(ucecSignatures['betterId']))]
mutInfoDfUCECMuts = maf_analysis_utils.summarize_per_case_mutation_info_for_mafs(ucecMuts)

mutSimInfoDf = mutation_modeling_util.compare_observed_to_expected_mut_burden_data(mutInfoDfUCECMuts, mutFracDictVals)

meltedDf = pd.melt(mutSimInfoDf, id_vars = ['Tumor_Sample_Barcode', 'Nmut'], value_vars=['ratioObservedToExpected',
                   'nOncogenicSnps',
                   'simulatedNOncogenicSNPs'])
meltedDf['ratio'] = meltedDf.apply(lambda row: row['value'] if row['variable'] == 'ratioObservedToExpected' else None, axis=1)

    
meltedDf.to_csv('~/Desktop/dataForLocalPlotting/tcgaUCECObservedVsExpectedRatios.tsv', sep='\t', index=False)





