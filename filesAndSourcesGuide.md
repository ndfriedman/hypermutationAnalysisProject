
# Figure 1 
<br>
<br>
## Essential files for figure 1 

all files copied from /juno/work/ccs/resources/impact/cbio_mutations on nov 19th 2019
**MAF with oncogkb and hotspot annotations**: 
/juno/work/taylorlab/friedman/myAdjustedDataFiles/data_mutations_extended_annotated_nov19_2019.maf
/juno/work/taylorlab/friedman/myAdjustedDataFiles/data_mutations_extended_annotated_sigContext_nov19_2019.maf

add trinuc and quadnuc information with python /ifs/work/taylorlab/friedman/myUtils/mutationSigUtils.py --mode trinucOnly --inputMaf /ifs/work/taylorlab/friedman/myAdjustedDataFiles/impactMafs/data_mutations_unfiltered_reviewed_oncokb.txt --outputDir /ifs/work/taylorlab/friedman/myAdjustedDataFiles --outputFilename annotatedOncoPlusHotspotMafAllImpact_trinuc
followed by the mutationSigUtils.create_reference_four_nuc() function
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

All main figure plots are made with 'plotFigure1.R' and all supplementary plots are made with 'plotSupplementaryFiguresFig1.R'

#### figure 1a
Define the hypermutation thresholds with: hypermutationAnalysisProject/plottingScripts/plotAndDefineHypermutationThresholds.R
Adjust them with /Users/friedman/Desktop/WORK/hypermutationProjectJupyterScripts/scriptsToGenerateFigures/generate_mut_classification_figure#1a#.ipynb
FINAL OUTPUTS ARE SAVED IN: 
/juno/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/hypermutationStatusIds

#### figure 1b
Summarize cancer types and signatures with #generate_cancer_types_and_sigs_figure#figure1b_c

#### figure 1c
Summarize cancer types and signatures with #generate_cancer_types_and_sigs_figure#figure1b_c

#### figure 1d
Summarize n hotposts/oncogenic info with generate_n_oncogenic_by_hyper_not_hyper_#figure1d#

#### figure 1e
Get n expected rates for all cases with:
Get all possible mutation summaries with /juno/work/taylorlab/friedman/hypermutationAnalysisProj/mutSimulation/summarize_possible_mutations.py, which writes a file to: /juno/work/taylorlab/friedman/myAdjustedDataFiles/allPossibleIMPACTMutationsSummary.tsv

use unadjustedSignatures from: '/juno/work/taylorlab/friedman/myAdjustedDataFiles/impactSignatureCalls_Nov20_2019_not_merged.tsv'
and run this script (it takes a long time):
python /juno/work/taylorlab/friedman/hypermutationAnalysisProjmutSimulation/summarize_expected_data.py
data is written here: /juno/work/taylorlab/friedman/hypermutationAnalysisProj/mutSimulation/expectedMutationTables/allHypermutatorsExpectedGeneMutInfo.tsv

-----------------------------------------------------------------------------------------------------------------

# Figure 2

#### figure 2a
compare n mutations in cancer type related and unrelated genes (note we need to define what is and what is not related)

#### figure 2b
?gene by gene hypermutated vs normal

#### figure 2c
unfilteredMaf is here: /juno/work/taylorlab/friedman/myAdjustedDataFiles/dataMutationsUnfilteredNov19.txt
use make_dnds_figure_2c to write dnds files to: /juno/work/taylorlab/friedman/myAdjustedDataFiles/mafsForDndsCV/
Run DNDS using:
runDnDsCv.R

#### figure 2d

#### figure 2e

# Figure 3

# Figure 4

CLONALITY MAF:
get CNCF filename info with: 'python myUtils/create_cncf_or_rdata_file_list.py'
file is by default written to: /juno/work/taylorlab/friedman/myAdjustedDataFiles/cncf_filenames.txt

Using: myUtils/runAnnotateMaf.R ,file is written to: /juno/work/taylorlab/friedman/myAdjustedDataFiles/filteredMaf_Nov19_2019_withCNCFAnnotation.maf


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



