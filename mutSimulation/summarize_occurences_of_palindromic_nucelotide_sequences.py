#written by noah friedman
#the goal of this code is to enumerate the places in the genome where:
#the nucleotide sequence TCGA occurs where the CGA codes for Arg
#this is because the R* and R->Q mutations are equally likely to occur in different directions
#so it makes it a fruitful base for analysis

import os
import pandas as pd
from collections import Counter
import re


tcgaPentaNucleotides = ['ATCGA', 'GTCGA', 'CTCGA', 'TTCGA',
						'TCGAA', 'TCGAC', 'TCGAG', 'TCGAT']

cntr = 0
mainDir = '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/simulatedMafs/allPossibleGeneMutMafWithPentanucleotideContext'
listOfDataFrames = []
for file in os.listdir(mainDir):
	localDict = {}
	gene = file.split('_')[0]
	
	#geneDf = pd.read_table(os.path.join(mainDir, file))
	geneDf = pd.read_csv(os.path.join(mainDir, file),
        header=0,
        usecols=["Ref_Tri.1", "Hugo_Symbol", "Consequence", "is-a-hotspot", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Type", "HGVSp_Short"], sep='\t')

	geneDfPenta = geneDf[geneDf['Ref_Tri.1'].isin(tcgaPentaNucleotides)]
	geneDfPenta['refAminoAcid'] = geneDfPenta['HGVSp_Short'].apply(lambda x:
		x[2] if isinstance(x, basestring) else None)
	geneDfPenta['altAminoAcid'] = geneDfPenta['HGVSp_Short'].apply(lambda x:
		x[len(x) - 1] if isinstance(x, basestring) else None)
	geneDfPenta['aminoAcidPosition'] = geneDfPenta['HGVSp_Short'].apply(lambda x:
		re.findall(r'\d+', x)[0] if isinstance(x, basestring) else None)
	geneDfRefAndContext = geneDfPenta[(geneDfPenta['refAminoAcid'] == 'R') & ((geneDfPenta['altAminoAcid'] == 'Q') | (geneDfPenta['altAminoAcid'] == '*'))]
	listOfDataFrames.append(geneDfRefAndContext)

	print cntr, gene
	cntr += 1

fullDf = pd.concat(listOfDataFrames)

fullDf.to_csv('/ifs/work/taylorlab/friedman/myAdjustedDataFiles/palindromicPoleContextSummary.tsv', index=False, sep='\t')
