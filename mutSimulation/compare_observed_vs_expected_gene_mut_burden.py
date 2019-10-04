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

msiProbByGene = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/relativeMsiProbabilityByGene.tsv')

geneDistributionsDf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/impact_gene_trinuc_distributions.tsv')


#SET up per gene snp rates
impactMaf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/tempScriptFiles/filteredMafAnnotated_trinuc.maf')
impactMaf['quadNuc'] = impactMaf.apply(lambda row: mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)

possibleMutData = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/simulatedMafs/allmut_impact_simulated_muts_mafAnno_trinuc.maf')
oncogenicSusceptibilityDict = mutation_modeling_util.calculate_quadnuc_based_oncogenic_susceptibility_dict(possibleMutData)

mutation_modeling_util.calculate_quadnuc_and_gene_based_oncogenic_susceptibility_dict(oncogenicSusceptibilityDict, possibleMutData)


Counter(possibleMutData['Hugo_Symbol'])
for hugo in set(possibleMutData['Hugo_Symbol']):
    print possibleMutData[possibleMutData['Hugo_Symbol'] == hugo].shape

