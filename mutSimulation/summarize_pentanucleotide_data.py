#wrriten by Noah Friedman
import os
import pandas as pd
from collections import Counter

def create_strand_specific_pentanucleotide_change(refPenta, refAllele, altAllele, variantType):

	complementDict = {'C': 'G', 'G': 'C', 'T': 'A', 'A': 'T'}

	if variantType != 'SNP': return None
	if not isinstance(refPenta, basestring): return None
	if len(refPenta) < 5: return None

	if refAllele != refPenta[2]:
		altAllele = complementDict[altAllele]
	#if refAllele != refPenta[2]: #this means the ref penta is corrected
		#refPenta = refPenta[::-1] #reverse the reference pentanucleotide so that it matches
		#refPenta = ''.join([complementDict[x] for x in refPenta])#and complement it
	return refPenta[:2] + '(' + refPenta[2] + '>' + altAllele + ')' + refPenta[3:]



cntr = 0
mainDir = '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/simulatedMafs/allPossibleGeneMutMafWithPentanucleotideContext'
listOfDicts = []
for file in os.listdir(mainDir):
	localDict = {}
	gene = file.split('_')[0]
	
	#geneDf = pd.read_table(os.path.join(mainDir, file))
	geneDf = pd.read_csv(os.path.join(mainDir, file),
        header=0,
        usecols=["Ref_Tri.1", "Hugo_Symbol", "Consequence", "is-a-hotspot", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Type"], sep='\t')

	if geneDf.shape[0] > 1:

		geneDf['pentaChange'] = geneDf.apply(lambda row: 
			create_strand_specific_pentanucleotide_change(row['Ref_Tri.1'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)

		nonSilentDf = geneDf[~geneDf['Consequence'].isin(set(['synonymous_variant', 'splice_region_variant,synonymous_variant', 'intron_variant', 'splice_region_variant,intron_variant']))]
		truncatingDf = geneDf[geneDf['Consequence'] == 'stop_gained']
		hotspotDf = geneDf[geneDf['is-a-hotspot'].notnull()]


		nonSilentPentaMutsDict = dict(nonSilentDf['pentaChange'].value_counts())
		truncatingPentaMutsDict = dict(truncatingDf['pentaChange'].value_counts())
		hotspotPentaMutsDict = dict(hotspotDf['pentaChange'].value_counts())

		localDict['Hugo_Symbol'] = gene
		localDict['nNonSilent'] = sum(nonSilentPentaMutsDict.values())
		localDict['nTruncating'] = sum(truncatingPentaMutsDict.values())
		localDict['nHotspot'] = sum(hotspotPentaMutsDict.values())

		for key,value in nonSilentPentaMutsDict.items():
			localDict[key+'_nonSilent'] = value

		for key,value in truncatingPentaMutsDict.items():
			localDict[key+'_truncating'] = value

		for key,value in hotspotPentaMutsDict.items():
			localDict[key+'_hotspot'] = value

		listOfDicts.append(localDict)

	print cntr, gene
	cntr += 1

df = pd.DataFrame(listOfDicts)

df.to_csv('/ifs/work/taylorlab/friedman/myAdjustedDataFiles/pentaNucMutationSummary.tsv', index=False, sep='\t')









