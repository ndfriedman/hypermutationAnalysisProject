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


    
    
    
def model_oncogenic_mut_acquisition(allSpectraDict, qNucDict, qNucDictAdj, geneOncogenicOrHotspotDict, geneMutProbInfo, simulatedNmuts = 1000, simulatedKCases = 20, mutationMode='oncogenic'):
    
    #initiate the signatures dicts
    signatureCols = ['Signature.' + str(i) for i in range(1,31)]
    sigsDict = dict()
    for col in signatureCols:
        sigsDict[col] = 0
     
    simulatedNmuts = 1000
    simulatedKCases = 20
    simulatedVals = []
    for col in signatureCols:
        sigsDict[col] = 1
        val1 = mutation_modeling_util.do_k_simulations_of_mutations(simulatedKCases,                                      
                                      sigsDict, 
                                      allSpectraDict, #a dictionary of all spectras from which we get a mutation signature spectra
                                      simulatedNmuts, #the number of mutations we are supposed to simulate
                                      qNucDict, #probability of picking a specific quadnuc given a spectra 
                                      qNucDictAdj,
                                      geneOncogenicOrHotspotDict, #probability of an oncogenic mutation O given a mutation in a gene G and a quadnuc Q
                                      geneMutProbInfo, #The probability of mutating a gene (based on the probabilities of genes being mutated across all silent mutations)
                                      mutMode=mutationMode, #mutMode--do we care about picking oncogenic mutations or picking hotspot mutations
                                      modelName='model1' #model: specifies the assumptions we use to pick the mutations
                                      )
        
        val2 = mutation_modeling_util.do_k_simulations_of_mutations(simulatedKCases,                                      
                                      sigsDict, 
                                      allSpectraDict, #a dictionary of all spectras from which we get a mutation signature spectra
                                      simulatedNmuts, #the number of mutations we are supposed to simulate
                                      qNucDict, #probability of picking a specific quadnuc given a spectra 
                                      qNucDictAdj,
                                      geneOncogenicOrHotspotDict, #probability of an oncogenic mutation O given a mutation in a gene G and a quadnuc Q
                                      geneMutProbInfo, #The probability of mutating a gene (based on the probabilities of genes being mutated across all silent mutations)
                                      mutMode=mutationMode, #mutMode--do we care about picking oncogenic mutations or picking hotspot mutations
                                      modelName='model2' #model: specifies the assumptions we use to pick the mutations
                                      )
        
        val3 = mutation_modeling_util.do_k_simulations_of_mutations(simulatedKCases,                                      
                                      sigsDict, 
                                      allSpectraDict, #a dictionary of all spectras from which we get a mutation signature spectra
                                      simulatedNmuts, #the number of mutations we are supposed to simulate
                                      qNucDict, #probability of picking a specific quadnuc given a spectra 
                                      qNucDictAdj,
                                      geneOncogenicOrHotspotDict, #probability of an oncogenic mutation O given a mutation in a gene G and a quadnuc Q
                                      geneMutProbInfo, #The probability of mutating a gene (based on the probabilities of genes being mutated across all silent mutations)
                                      mutMode=mutationMode, #mutMode--do we care about picking oncogenic mutations or picking hotspot mutations
                                      modelName='model3' #model: specifies the assumptions we use to pick the mutations
                                      )
        
        simulatedVals.append((col, val1, val2, val3))
        sigsDict[col] = 0
        
    return simulatedVals

    

mafHere = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/annotatedOncoPlusHotspotMafAllImpact_trinuc')
#mafHere = pd.read_table('~/Desktop/OfflineStuffFor11_20Trip/annotatedOncoPlusHotspotMafAllImpact_trinuc')
#geneDist = pd.read_table('~/Desktop/OfflineStuffFor11_20Trip/impact_gene_trinuc_distributions.tsv')

mafHere['quadNuc'] = mafHere.apply(lambda row: mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)

#the background distribtuion of all observed mutations in impact
mafBackground = pd.read_table(pathPrefix + '/ifs/res/taylorlab/ang46/ext/dmp/mskimpact/mutation_data.txt')

#The background distribution of trinucs 
geneDistributionsDf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/impact_gene_trinuc_distributions.tsv')

#smallMafWithOncogenicAnnotations = mafHere.sample(3000, axis=0)
#smallMafBackgroundMaf = mafBackground.sample(3000, axis=0)
#smallMafWithOncogenicAnnotations['quadNuc'] = smallMafWithOncogenicAnnotations['quadNuc'] = smallMafWithOncogenicAnnotations.apply(lambda row: 
#   			mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)

    
#BIG LINE OF CODE THAT GETS ALL THE PROBABILITY MODELS WE NEED
#DO it quick with minimafs
#quadNucDict, quadNucDictAdjusted, geneOncogenicOrHotspotProbDict, geneMutProbInfo = mutation_modeling_util.initiate_models(smallMafWithOncogenicAnnotations, smallMafBackgroundMaf, geneDistributionsDf)
#Do it the long way with big ass full mafs
quadNucDict, quadNucDictAdjusted, geneOncogenicOrHotspotProbDict, geneMutProbInfo = mutation_modeling_util.initiate_models(mafHere, mafBackground, geneDistributionsDf)
mutation_modeling_util.write_model_information(quadNucDict=quadNucDict, quadNucDictAdjusted=quadNucDictAdjusted, geneOncogenicOrHotspotProbDict=geneOncogenicOrHotspotProbDict, geneMutProbInfo=geneMutProbInfo, write_dir='~/Desktop/mnt/ifs/work/taylorlab/friedman/myAdjustedDataFiles/mutationSimultationData',)

spectraD = mutationSigUtils.convert_spectrum_file_to_dict_of_dicts(spectrumFile=pathPrefix + '/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt')

#DO SIMULATION
vals = model_oncogenic_mut_acquisition(spectraD, quadNucDict, quadNucDictAdjusted, geneOncogenicOrHotspotProbDict, geneMutProbInfo, simulatedNmuts = 1000, simulatedKCases = 50, mutationMode='hotspot')

listOfDicts = []
#MAKE a summary dataframe using the long data format
for val in vals:
    localDict = dict()
    localDict['SigName'] = val[0]
    localDict['average'] = val[1]
    localDict['modelName'] = 'model1'
    listOfDicts.append(localDict)
    
    localDict = dict()
    localDict['SigName'] = val[0]
    localDict['average'] = val[2]
    localDict['modelName'] = 'model2'
    listOfDicts.append(localDict)
    
    localDict = dict()
    localDict['SigName'] = val[0]
    localDict['average'] = val[3]
    localDict['modelName'] = 'model3'
    listOfDicts.append(localDict)

dfData = pd.DataFrame(listOfDicts)

dfData.to_csv('~/Desktop/dataForLocalPlotting/simulationData.tsv', sep='\t', index=False)






quadNucDict, geneOncogenicOrHotspotProbDict = mutation_modeling_util.do_initiation(mafWithInfo, geneDistributions=pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/impact_gene_trinuc_distributions.tsv')


l = []
for i in range(100):
    l.append(mutation_modeling_util.pick_mutation_given_gene_and_quadNuc('ARID1A', 'TTCG', geneOncogenicOrHotspotProbDict, mode='hotspot'))
print Counter(l)


tmzSigs = gliomaSigs[gliomaSigs['mean_11'] > .25]
renameDict = {key:value for (key,value) in [('mean_' + str(i), 'Signature.' + str(i)) for i in range(1,31)]}
tmzSigs = tmzSigs.rename(columns=renameDict)
signatureCols = ['Signature.' + str(i) for i in range(1,31)]


val = asses_simulated_mutation_rates_by_signatures(sigMagnitudesDict, 100, spectraD, 1000, quadNucDict, geneOncogenicOrHotspotProbDict, mMode='oncogenic')
df = pd.DataFrame([])
df['signature'] = val.keys()
df['average'] = val.values()

df.to_csv('~/Desktop/noahQuickTest.tsv', sep='\t', index=False)

cases = set(tmzSigs['Tumor_Sample_Barcode'])
for case in cases:
    print mutDataDf[mutDataDf['Tumor_Sample_Barcode'] == case]['nMuts'].iloc[0], tmzSigs[tmzSigs['Tumor_Sample_Barcode'] == case]['Nmut'].iloc[0]
    #localD = tmzSigs[tmzSigs['Tumor_Sample_Barcode'] == case]
    #sigMagnitudesDict = localD[signatureCols].iloc[0].to_dict()
    
    #sigMagnitudesDict = {key:value for (key, value) in [('Signature.' + str(i), 0) for i in range(1,31)]}
    #sigMagnitudesDict['Signature.1'] = 1
    
    #mutation_modeling_util.do_k_simulations_of_mutations(sigMagnitudesDict, 10, spectraD, 1000, quadNucDict, geneOncogenicOrHotspotProbDict, mutMode='oncogenic')
    #spectraChosen, quadNucsChosen, genesChosen, mutationsChosen = mutation_modeling_util.pick_n_mutations_given_spectra(sigMagnitudesDict, spectraD, 1000, quadNucDict, geneOncogenicOrHotspotProbDict, mutMode='oncogenic')
    #print Counter(spectraChosen), Counter(quadNucsChosen), Counter(genesChosen), Counter(mutationsChosen)
    #break
    
    
#################################################################
#########Look at mutation landscapes across different mutation signatures
    
dfTest = pd.read_table(pathPrefix + '/ifs/res/taylorlab/ang46/ext/dmp/mskimpact/mutation_data.txt')


