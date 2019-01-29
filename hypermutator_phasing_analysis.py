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

phasingData = pd.read_table(pathPrefix + '/home/ang46/luna/projects/collaborations/noah/phased_samples_for_noah_20190116.txt')


phasingData.columns.values
phasingData['oncogenic.2']

oneOncMutant = phasingData[(phasingData['oncogenic.1'].notnull()) | (phasingData['oncogenic.2'].notnull())]

oneOncMutant.columns.values
tenMostCommonMuts = Counter(oneOncMutant['Hugo_Symbol']).most_common(10)
for mut in tenMostCommonMuts:
    gene = mut[0]
    geneMuts = oneOncMutant[oneOncMutant['Hugo_Symbol'] == gene]
    print gene, geneMuts[geneMuts['phase'] != 'unknown']['phase']

