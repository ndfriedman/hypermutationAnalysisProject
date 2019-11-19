MAF with oncogkb and hotspot annotations: 
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



