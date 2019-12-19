#written by noah friedman
import os
import pandas as pd
from collections import Counter
import sys
sys.path.append('/ifs/work/taylorlab/friedman/myUtils')
import mutationSigUtils 

mode = sys.argv[1]
print 'mode is ', mode
cntr = 0
mainDir = '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/simulatedMafs/allPossibleGeneMutMafWithPentanucleotideContext'
listOfDfs = []
for file in os.listdir(mainDir):
	localDict = {}
	gene = file.split('_')[0]
	print gene, cntr
	geneDf = pd.read_csv(os.path.join(mainDir, file),
        header=0,
        usecols=["Ref_Tri.1", "Hugo_Symbol", "Consequence", "is-a-hotspot", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Type", "Ref_Tri", "HGVSp_Short", "Start_Position"], sep='\t')


	if geneDf.shape[0] > 1:

		if mode == 'hotspot':
			gdf = geneDf[geneDf['is-a-hotspot'].notnull()]
			gdf = gdf[gdf['is-a-hotspot'] == 'Y']
		elif mode == 'truncating':
			gdf = geneDf[geneDf['Consequence'] == 'stop_gained']

		if gdf.shape[0] > 1:

			gdf['quadNuc'] = gdf.apply(lambda row: mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)
			gdf['allele'] = gdf['Hugo_Symbol'] + '_' + gdf['HGVSp_Short']
			
			listOfDfs.append(gdf)
	cntr += 1

fullDf = pd.concat(listOfDfs)
print 'wriiting to file'

outfile = ''
if mode == 'hotspot':
	outfile = '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/hotspotContextSummary.tsv'
if mode == 'truncating':
	outfile = '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/truncatingContextSummary.tsv'

fullDf.to_csv(outfile, index=False, sep='\t')
