
# Figure 1 
<br>
<br>
## Essential files for figure 1 

all files copied from /juno/work/ccs/resources/impact/cbio_mutations on nov 19th 2019
**MAF with oncogkb and hotspot annotations**: 
/juno/work/taylorlab/friedman/myAdjustedDataFiles/data_mutations_extended_annotated_nov19_2019.maf
37990 cases; 331821 lines; 

**Text file with cancer type information**
/juno/work/taylorlab/friedman/myAdjustedDataFiles/cancerTypeInfo_asOfNov192019.txt
37989 cases, access with analysis_utils.get_cancer_type_information()

**Text file with TMB and MSI stats**
/juno/work/taylorlab/friedman/myAdjustedData/mutations_TMB_and_MSI_stats.txt
37990 cases

**Text file with IMPACT signature decompositions from Alex Penson**
/juno/work/taylorlab/friedman/myAdjustedDataFiles/impactSignatureCalls_Nov20_2019.tsv
--note the tmb and nmut values alex generated were overwritten and replaced with stuff from chai
--I added information about dominant signatures here too

#### STEP 1: figure 1a
Define the hypermutation thresholds with: hypermutationAnalysisProject/plottingScripts/plotAndDefineHypermutationThresholds.R
Adjust them with /Users/friedman/Desktop/WORK/hypermutationProjectJupyterScripts/scriptsToGenerateFigures/generate_mut_classification_figure#1a#.ipynb
FINAL OUTPUTS ARE SAVED IN: 
/juno/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/hypermutationStatusIds

#### STEP 2: figure 1b/c
Summarize cancer types and signatures with 


<br>
<br>
<br>
<br>
<br>
<br>
#OLD

/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/all_impact_mutations_annotated_cohort.maf
275,171 lines
 
MAFs with facets annotations from MAF anno:
/ifs/work/taylorlab/friedman/myAdjustedDataFiles/subsettedMafs/Endometrial_HypermutantCaseMuts_MAF_ANNO.maf
/ifs/work/taylorlab/friedman/myAdjustedDataFiles/subsettedMafs/Colorectal_HypermutantCaseMuts_MAF_ANNO.maf
/ifs/work/taylorlab/friedman/myAdjustedDataFiles/subsettedMafs/Glioma_HypermutantCaseMuts_MAF_ANNO.maf


SIGNATURES
Signature 30 spectrums: /ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt
Signatures from Alex Penson: /ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/signatures_from_unfiltered_maf.txt
TCGA signatures: /home/friedman/junoData/myAdjustedDataFiles/tcgaSigsCombined.txt


HYPERMUTATION CLASSIFICATIONS:
/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/hypermutationStatusIds

SIMULATED ALL MUTATION MAFS:
with pentanucleotide context: /ifs/work/taylorlab/friedman/myAdjustedDataFiles/simulatedMafs/allPossibleGeneMutMafWithPentanucleotideContext
summary:
/ifs/work/taylorlab/friedman/myAdjustedDataFiles/pentaNucMutationSummary.tsv

PENTA nucleotide summaries:
/ifs/work/taylorlab/friedman/myAdjustedDataFiles/pentaNucMutationSummary.tsv

ALL MUTS SUMMARY:
use python summarize_possible_mutations.py to generate:
/juno/work/taylorlab/friedman/myAdjustedDataFiles/allPossibleIMPACTMutationsSummary.tsv



EXPECTED RATES
uses /juno/work/taylorlab/friedman/myAdjustedDataFiles/allPossibleIMPACTMutationsSummary.tsv as an input for all numbers of mutations, uses the function: mutation_modeling_util.convert_counts_information_to_fraction(countsDf) to get intermediate data frame of chances for mutations
expected truncating, hotspot and oncogenic rates by gene for all hypermutated cases /ifs/work/taylorlab/friedman/hypermutationAnalysisProj/mutSimulation/expectedMutationTables/allHypermutatorsExpectedGeneMutInfo.tsv
expected rates for tcga hypermutated cases: /ifs/work/taylorlab/friedman/hypermutationAnalysisProj/mutSimulation/expectedMutationTables/tcgaHypermutatorsExpectedGeneMutInfo.tsv

TCGA





ESOTERIC:
palindromic pole contexts: /ifs/work/taylorlab/friedman/myAdjustedDataFiles/palindromicPoleContextSummary.tsv



