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
import math


def get_average_simulated_n_oncogenic_muts_for_case(
        simulatedKCases, 
        sigsDict,
        allSpectraDict,
        simulatedNmuts,
        qNucDict,
        qNucDictAdj,
        geneOncogenicOrHotspotDict,
        geneMutProbInfo
        ):
    
    mutationMode = 'oncogenic'
    
    val1, simMuts1 = mutation_modeling_util.do_k_simulations_of_mutations(simulatedKCases,                                      
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
        
    val2, simMuts2 = mutation_modeling_util.do_k_simulations_of_mutations(simulatedKCases,                                      
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
        
    val3, simMuts3 = mutation_modeling_util.do_k_simulations_of_mutations(simulatedKCases,                                      
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
        
    return np.nanmean([val1, val2, val3]), [simMuts1, simMuts2, simMuts3]
    

def simulate_hotspot_acquisition_across_tmz_cases(tmzCaseSigs,
                                                  qNucDict,
                                                  qNucDictAdj,
                                                  gOncogenicOrHotspotProbDict,
                                                  gMutProbInfo,
                                                  spectraD
                                                  ):
    
    simulatedOncogenicMuts = {}
    signatureCols = ['Signature.' + str(i) for i in range(1,31)]
    
    for index, row in tmzCaseSigs.iterrows():
        simulatedNmuts = 1000
        simulatedKCases = 20
        
        case = row['Tumor_Sample_Barcode']
        print 'doing mutation simulation for case: ', case
        curCaseSigs = row[signatureCols].to_dict()
        
        v, simulatedMuts = get_average_simulated_n_oncogenic_muts_for_case(
                simulatedKCases, 
                curCaseSigs,
                spectraD,
                simulatedNmuts,
                qNucDict,
                qNucDictAdj,
                gOncogenicOrHotspotProbDict,
                gMutProbInfo
        )
        
        simulatedOncogenicMuts[case] = [v, simulatedMuts]
    return simulatedOncogenicMuts



def compare_observed_vs_simulated_mut_burden_tmz(nSnpsD, nOncogenicSnpsD, simD):
    for key, value in simD.items():
        #print key, value
        nSnpsInCase = nSnpsD[key]
        nOncogenicSnpsInCase = nOncogenicSnpsD[key]
        nSimOncogenicSnpsForCase = value*nSnpsInCase
        print nOncogenicSnpsInCase/nSimOncogenicSnpsForCase

def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z


def get_simulated_mutation_summary(simData, normalize=True):
    
    def summarize_n_times_gene_mutated_per_case(data):
        genes = []
        for entry in data:
            for mutation in entry:
                gene = mutation.split('_')[0]
                genes.append(gene)
        
        #have the chance that any one mutation = any one gene
        return {key: value*1.0/(len(data)*1000 )for (key, value) in Counter(genes).items()}
    
    #if calls simplifies the simFracDict so that each Tumor Sample Barcode maps to 1 average times a gene is mutated per mutation
    def normalize_sim_frac_dict(sfd):
        impactGenes = set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KMT2B', 'KMT5A', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51C', 'RAD51L1', 'RAD51L3', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2'])
        normedDict = {}
        for key, value in sfd.items():
            vals = []
            localDict = {}
            for gene in impactGenes:
                if gene in value[0]:
                    vals.append(value[0][gene])
                else:
                    vals.append(0)
                
                if gene in value[1]:
                    vals.append(value[1][gene])
                else:
                    vals.append(0)
                    
                if gene in value[2]:
                    vals.append(value[2][gene])
                else:
                    vals.append(0)
                    
                localDict[gene] = np.nanmean(vals)
            normedDict[key] = localDict
        return normedDict
    
    
    
    simFracDict = {}
    for key, value in simData.items():
        model1Data = value[1][0]
        model2Data = value[1][1]
        model3Data = value[1][2]
        
        perCaseMutInfo1 = summarize_n_times_gene_mutated_per_case(model1Data)
        perCaseMutInfo2 = summarize_n_times_gene_mutated_per_case(model2Data)
        perCaseMutInfo3 = summarize_n_times_gene_mutated_per_case(model3Data)
        simFracDict[key] = [perCaseMutInfo1, perCaseMutInfo2, perCaseMutInfo3]

    return normalize_sim_frac_dict(simFracDict)

def get_n_simulated_gene_muts_per_case(simData, dfWithNmut):
    impactGenes = set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KMT2B', 'KMT5A', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51C', 'RAD51L1', 'RAD51L3', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2'])
    geneSimDict = dict()
    for gene in impactGenes:
        geneSimDict[gene] = 0
    
    for key, value in simData.items():
        nmutForCase = dfWithNmut[dfWithNmut['Tumor_Sample_Barcode'] == key]['Nmut'].iloc[0]
        for gene, ratio in value.items():
            geneSimDict[gene] = geneSimDict[gene] + nmutForCase*ratio

    for key, value in geneSimDict.items():
        geneSimDict[key] = value/len(simData)

    return geneSimDict

###############SETUP########################

quadNucDict =  mutation_modeling_util.do_initiation(geneDistributions=pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/impact_gene_trinuc_distributions.tsv')
spectraD = mutationSigUtils.convert_spectrum_file_to_dict_of_dicts(spectrumFile=pathPrefix + '/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt')

mafHere = pd.read_table('~/Desktop/OfflineStuffFor11_20Trip/annotatedOncoPlusHotspotMafAllImpact_trinuc')
mafBackground = pd.read_table(pathPrefix + '/ifs/res/taylorlab/ang46/ext/dmp/mskimpact/mutation_data.txt')
geneDistributionsDf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/impact_gene_trinuc_distributions.tsv')
mafHere['quadNuc'] = mafHere.apply(lambda row: 
   	mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)

#Do it the long way with big ass full mafs
quadNucDict, quadNucDictAdjusted, geneOncogenicOrHotspotProbDict, geneMutProbInfo = mutation_modeling_util.initiate_models(mafHere, mafBackground, geneDistributionsDf)

############################################

impactSigs = pd.read_table(pathPrefix + '/ifs/res/taylorlab/impact_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
impactSigs['pid'] = impactSigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
impactSigs['cancer_type'] = impactSigs['pid'].apply(lambda x: cDict[x] if x in cDict else None)
gliomaSigs = impactSigs[impactSigs['cancer_type'] == 'Glioma']

tmzSigs = gliomaSigs[(gliomaSigs['mean_11'] > .25) & (gliomaSigs['Nmut'] > 50)]
nonTMZGliomaSigs = gliomaSigs[(gliomaSigs['mean_11'] <= .25) & (gliomaSigs['Nmut'] <= 50)]
renameDict = {key:value for (key,value) in [('mean_' + str(i), 'Signature.' + str(i)) for i in range(1,31)]}
tmzSigs = tmzSigs.rename(columns=renameDict)
nonTMZGliomaSigs = nonTMZGliomaSigs.rename(columns=renameDict)

caseSimDict = simulate_hotspot_acquisition_across_tmz_cases(tmzSigs,
                                              quadNucDict,
                                              quadNucDictAdjusted,
                                              geneOncogenicOrHotspotProbDict,
                                              geneMutProbInfo,
                                              spectraD
                                              )


caseSimDictNonTMZ = simulate_hotspot_acquisition_across_tmz_cases(nonTMZGliomaSigs,
                                              quadNucDict,
                                              quadNucDictAdjusted,
                                              geneOncogenicOrHotspotProbDict,
                                              geneMutProbInfo,
                                              spectraD
                                              )

##############################################


#TODO compare simulations graph

tmzCasesMuts = mafHere[mafHere['Tumor_Sample_Barcode'].isin(set(tmzSigs['Tumor_Sample_Barcode']))]
nonTMZCasesMuts = mafHere[mafHere['Tumor_Sample_Barcode'].isin(set(nonTMZGliomaSigs['Tumor_Sample_Barcode']))]

gliomaMuts = mafHere[mafHere['Tumor_Sample_Barcode'].isin(set(gliomaSigs['Tumor_Sample_Barcode']))]
gliomaMuts['pid'] = gliomaMuts['Tumor_Sample_Barcode'].apply(lambda x: x[:9])

allCaseMuts = mafHere[mafHere['Tumor_Sample_Barcode'].isin(set(nonTMZGliomaSigs['Tumor_Sample_Barcode']) | set(tmzSigs['Tumor_Sample_Barcode']))]
mutInfoDf = maf_analysis_utils.summarize_per_case_mutation_info_for_mafs(allCaseMuts)
mutInfoDf['isTMZ'] = mutInfoDf['Tumor_Sample_Barcode'].apply(lambda x: True if x in caseSimDict else False)

mergedDict = merge_two_dicts(caseSimDict, caseSimDictNonTMZ)

mutInfoDf['Nmut']

mutInfoDf['nSimulatedOncogenicMutRatio'] = mutInfoDf['Tumor_Sample_Barcode'].apply(lambda x: mergedDict[x])
mutInfoDf['nSimulatedOncogenicMuts'] = mutInfoDf.apply(lambda row: row['Nmut'] * row['nSimulatedOncogenicMutRatio'], axis=1)
mutInfoDf['observedToSimulatedRatio'] = mutInfoDf.apply(lambda row: math.log(row['nOncogenicSnps']/row['nSimulatedOncogenicMuts'], 2) if row['nOncogenicSnps'] > 0 else None, axis = 1)

mutInfoDf= mutInfoDf.rename(columns={'observedToSimulatedRatio': 'Ratio of Observed Oncogenic muts/Simulated Oncogenic Muts',
                          'nSimulatedOncogenicMuts': 'Number of simulated Oncogenic muts',
                          'nOncogenicSnps': 'Number of Oncogenic Snps'})

meltedDf = pd.melt(mutInfoDf, id_vars = ['Tumor_Sample_Barcode', 'Nmut'], value_vars=['Ratio of Observed Oncogenic muts/Simulated Oncogenic Muts',
                   'Number of simulated Oncogenic muts',
                   'Number of Oncogenic Snps'])
    
meltedDf['ratio'] = meltedDf.apply(lambda row: row['value'] if row['variable'] == 'Ratio of Observed Oncogenic muts/Simulated Oncogenic Muts' else None, axis=1)

meltedDf.to_csv('~/Desktop/dataForLocalPlotting/simulatedVsObservedDataTMZ.tsv', sep='\t', index=False)


#######################################333
#Do stuff for comparing the number of mutations observed in a gene vs expected

caseSimDictWithMutSims = simulate_hotspot_acquisition_across_tmz_cases(tmzSigs,
                                              quadNucDict,
                                              quadNucDictAdjusted,
                                              geneOncogenicOrHotspotProbDict,
                                              geneMutProbInfo,
                                              spectraD
                                              )

sFD = get_simulated_mutation_summary(caseSimDictWithMutSims)
simPerCase = get_n_simulated_gene_muts_per_case(sFD, mutInfoDf)

occurenceCounter, rankingDict, gliomaGenesFractionalDict = maf_analysis_utils.enumerate_top_n_oncogenic_mutated_genes_across_cohort(gliomaMuts)

listOfDicts = []
for gene, value in gliomaGenesFractionalDict.items():
    if gene in simPerCase:
        listOfDicts.append({'Gene': gene, 'Ratio': math.log(value/simPerCase[gene], 2)})
    else:
        listOfDicts.append({'Gene': gene, 'Ratio': None})

df = pd.DataFrame(listOfDicts)
df.to_csv('~/Desktop/dataForLocalPlotting/tmzSimVsObservedGeneRatios.tsv', sep='\t', index=False)

############################################3
math.log(8, 2)

meltedDf['variable']

mutInfoDf.columns.values
#mafBackground = pd.read_table(pathPrefix + '/ifs/res/taylorlab/ang46/ext/dmp/mskimpact/mutation_data.txt')
#tmzCasesMuts = mafBackground[mafBackground['Tumor_Sample_Barcode'].isin(set(tmzSigs['Tumor_Sample_Barcode']))]

"""mutInfoDfTMZ = maf_analysis_utils.summarize_per_case_mutation_info_for_mafs(tmzCasesMuts)

nSnpDictTMZ = dict(zip(mutInfoDfTMZ['Tumor_Sample_Barcode'],mutInfoDfTMZ['NSnps']))
nOncogenicSnpsDictTMZ = dict(zip(mutInfoDfTMZ['Tumor_Sample_Barcode'], mutInfoDfTMZ['nOncogenicSnps']))
compare_observed_vs_simulated_mut_burden_tmz(nSnpDictTMZ, nOncogenicSnpsDictTMZ, caseSimDict)

mutInfoDfNonTMZ = maf_analysis_utils.summarize_per_case_mutation_info_for_mafs(nonTMZCasesMuts)
nSnpDictNonTMZ = dict(zip(mutInfoDfNonTMZ['Tumor_Sample_Barcode'],mutInfoDfNonTMZ['NSnps']))
nOncogenicSnpsDictNonTMZ = dict(zip(mutInfoDfNonTMZ['Tumor_Sample_Barcode'], mutInfoDfNonTMZ['nOncogenicSnps']))
compare_observed_vs_simulated_mut_burden_tmz(nSnpDictNonTMZ, nOncogenicSnpsDictNonTMZ, caseSimDictNonTMZ)
"""


