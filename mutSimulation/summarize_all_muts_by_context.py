#written by Noah Friedman

#written by noah friedman
import os
import pandas as pd
from collections import Counter
import sys
sys.path.append('/ifs/work/taylorlab/friedman/myUtils')
import mutationSigUtils 


cntr = 0
mainDir = '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/simulatedMafs/allPossibleGeneMutMafWithPentanucleotideContext'
listOfDfs = []
for file in os.listdir(mainDir):
	localDict = {}
	gene = file.split('_')[0]
	print gene, cntr
	geneDf = pd.read_csv(os.path.join(mainDir, file),
        header=0,
        usecols=["Ref_Tri.1", "Hugo_Symbol", "Consequence", "is-a-hotspot", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Type", "Ref_Tri", "HGVSp_Short", "Start_Position", "oncogenic"], sep='\t')

	if geneDf.shape[0] > 1:
		geneDf['quadNuc'] = geneDf.apply(lambda row: mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)
		geneDf['allele'] = geneDf['Hugo_Symbol'] + '_' + geneDf['HGVSp_Short']
			
		listOfDfs.append(geneDf)
	cntr += 1

fullDf = pd.concat(listOfDfs)

outfile = '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/allMutsContextSummary.tsv'

print 'writing file to', outfile
fullDf.to_csv(outfile, index=False, sep='\t')

