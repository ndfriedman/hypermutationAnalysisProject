#written by noah friedman, help from craig to run
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)
library(dndscv)
library(data.table)
library(ggrepel)

library(devtools)

fix_maf_cols <- function(maf){
  mutations <- maf[Variant_Type %in% c("SNP", "INS", "DEL"),.
                   (sampleID = Tumor_Sample_Barcode,
                     chr = Chromosome,
                     pos = Start_Position,
                     ref = Reference_Allele,
                     alt = Tumor_Seq_Allele2, 
                     gene= Hugo_Symbol,
                     hgvs =HGVSp_Short
                   )]
  return(mutations)
}

endometrialHyperMaf <- read.table('/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myAdjustedDataFiles/endometrialHypermutatedMutationsUnfiltered.maf', sep = '\t', header=TRUE)
endometrialHyperMaf <- data.table(endometrialHyperMaf)
endometrialHyperMaf <- fix_maf_cols(endometrialHyperMaf)

colorectalHyperMaf <- read.table('/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myAdjustedDataFiles/colorectalHypermutatedMutationsUnfiltered.maf', sep = '\t', header=TRUE)
colorectalHyperMaf <- data.table(colorectalHyperMaf)
colorectalHyperMaf <- fix_maf_cols(colorectalHyperMaf)

gliomaHyperMaf <- read.table('/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myAdjustedDataFiles/gliomaHypermutatedMutationsUnfiltered.maf', sep = '\t', header=TRUE)
gliomaHyperMaf <- data.table(gliomaHyperMaf)
gliomaHyperMaf <- fix_maf_cols(gliomaHyperMaf)

#NORMALS

gliomaNormalMaf <- read.table('/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myAdjustedDataFiles/gliomaNormalMutationsUnfiltered.maf', sep = '\t', header=TRUE)
gliomaNormalMaf <- data.table(gliomaNormalMaf)
gliomaNormalMaf <- fix_maf_cols(gliomaNormalMaf)

endometrialNormalMaf <- read.table('/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myAdjustedDataFiles/endometrialNormalMutationsUnfiltered.maf', ,sep = '\t', header=TRUE)
endometrialNormalMaf <- data.table(endometrialNormalMaf)
endometrialNormalMaf <- fix_maf_cols(endometrialNormalMaf)

colorectalNormalMaf <- read.table('/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myAdjustedDataFiles/colorectalNormalMutationsUnfiltered.maf', sep = '\t', header=TRUE)
colorectalNormalMaf <- data.table(colorectalNormalMaf)
colorectalNormalMaf <- fix_maf_cols(colorectalNormalMaf)

gene_list <- union(union(endometrialMaf$Hugo_Symbol,colorectalMaf$Hugo_Symbol),  gliomaMaf$Hugo_Symbol)
gene_list341 <- c('ABL1', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BAP1', 'BARD1', 'BBC3', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CARD11', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79B', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSF1R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EP300', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERG', 'ESR1', 'ETV1', 'ETV6', 'EZH2', 'FAM175A', 'FAM46C', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXP1', 'FUBP1', 'GATA1', 'GATA2', 'GATA3', 'GNA11', 'GNAQ', 'GNAS', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3B', 'HNF1A', 'HRAS', 'ICOSLG', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INPP4A', 'INPP4B', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAPK1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MITF', 'MLH1', 'KMT2A', 'KMT2D', 'KMT2C', 'MPL', 'MRE11A', 'MSH2', 'MSH6', 'MTOR', 'MUTYH', 'MYC', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOR1', 'NF1', 'NF2', 'NFE2L2', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTRK1', 'NTRK2', 'NTRK3', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLE', 'PPP2R1A', 'PRDM1', 'PRKAR1A', 'PTCH1', 'PTEN', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAC1', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RUNX1', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SETD2', 'SF3B1', 'SH2D1A', 'SHQ1', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SOCS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SRC', 'STAG2', 'STK11', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TBX3', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP63', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'VHL', 'VTCN1', 'WT1', 'XIAP', 'XPO1', 'YAP1', 'YES1')
gene_list410 <- c('ABL1', 'ACVR1', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'EP300', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1',  'KMT2A', 'KMT2D', 'KMT2C', 'MPL', 'MRE11A', 'MSH2', 'MSH6', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPM1D', 'PPP2R1A', 'PPP6C', 'PRDM1', 'PRKAR1A', 'PTCH1', 'PTEN', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAD21', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RUNX1', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHQ1', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SOCS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'VEGFA', 'VHL', 'VTCN1', 'WT1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2')
gene_list468 <- c('ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'AMER1', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'GTF2I', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1','KMT2A', 'KMT2D', 'KMT2C', 'KMT2B', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 
                  'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SETD8', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2')

#PLACATE DNDS cv by getting rid of genes it does not like
gene_list410 <- setdiff(gene_list410, c("MAP3K14", "MYCL1", "RYBP", "KMT5A", "PAK5", "FAM123B"))
gene_list468 <- setdiff(gene_list468, c("MAP3K14", "MYCL1", "RYBP", "KMT5A", "PAK5", "FAM123B"))

#RUN DNDS on everybody
dndsEndometrialHyper <- dndscv(endometrialHyperMaf, gene_list = gene_list468, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, use_indel_sites=F)
dndsColorectalHyper <- dndscv(colorectalHyperMaf, gene_list = gene_list468, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, use_indel_sites=F)
dndsGliomaHyper <- dndscv(gliomaHyperMaf, gene_list = gene_list468, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, use_indel_sites=F)

dndsEndometrialNormal <- dndscv(endometrialNormalMaf, gene_list = gene_list468, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, use_indel_sites=F)
dndsColorectalNormal <- dndscv(colorectalNormalMaf, gene_list = gene_list468, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, use_indel_sites=F)
dndsGliomaNormal <- dndscv(gliomaNormalMaf, gene_list = gene_list468, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, use_indel_sites=F)

sel_cvEdnoH <- dndsEndometrialHyper$sel_cv
sel_cvEdnoN <- dndsEndometrialNormal$sel_cv
mergedEndo <- merge(sel_cvEdnoN, sel_cvEdnoH, suffixes = c(".Normal",".Hypermutated"), by='gene_name')

sel_cvColoH <- dndsColorectalHyper$sel_cv
sel_cvColoN <- dndsColorectalNormal$sel_cv
mergedColo <- merge(sel_cvColoN, sel_cvColoH, suffixes = c(".Normal",".Hypermutated"), by='gene_name')

sel_cvGliomaH <- dndsGliomaHyper$sel_cv
sel_cvGliomaN <- dndsGliomaNormal$sel_cv
mergedGlioma <- merge(sel_cvGliomaN, sel_cvGliomaH, suffixes = c(".Normal",".Hypermutated"), by='gene_name')

emptyTheme <- theme(axis.line = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())

make_dnds_plot <- function(dndsData, title){
  capThresh <- 1e-10
  plotThresh <- .01
  minCoord = .9
  maxCoord = 2-log10(capThresh)
  
  ggplot()+
    geom_text_repel(data = dndsData[(dndsData$qNormal < plotThresh) | (dndsData$qHypermutated < plotThresh),],
                    aes(x=1- log10(qNormal), y=1-log10(qHypermutated), label=gene_name))+
    geom_point(data = dndsData[(dndsData$qNormal >= plotThresh) & (dndsData$qHypermutated >= plotThresh),],
               aes(x=1- log10(qNormal), y=1-log10(qHypermutated)))+

    scale_x_continuous(limits=c(minCoord, maxCoord), breaks=c(1,2,3,6,10))+
    scale_y_continuous(limits=c(minCoord, maxCoord), breaks=c(1,2,3,6,10))+
    
    #BELLS and whistles for the plot
    geom_segment(aes(x=minCoord, xend=maxCoord, y=1- log10(.01), yend= 1- log10(.01)), colour='black', linetype=2)+
    geom_segment(aes(x=1- log10(.01), xend=1- log10(.01), y=minCoord, yend= maxCoord), colour='black', linetype=2)+
    
    geom_segment(aes(x=minCoord, xend=maxCoord, y=1- log10(.1), yend= 1- log10(.1)), colour='black')+
    geom_segment(aes(x=1- log10(.1), xend=1- log10(.1), y=minCoord, yend= maxCoord), colour='black')+
    ylab('1 minus log(q) value in hypermutators')+
    xlab('1 minus log(q) value in non-hypermutators')+
    emptyTheme+
    ggtitle(title)
}

make_dnds_weak_drivers_plot <- function(dndsData, title){
  plotThresh <- .1
  minCoord = .9
  
  dndsData <- dndsData[(dndsData$qglobal_cv.Normal > .01) & (dndsData$qglobal_cv.Normal > .01),]
  ggplot()+
    geom_text_repel(data = dndsData[(dndsData$qglobal_cv.Normal < plotThresh) | (dndsData$qglobal_cv.Hypermutated < plotThresh),],
                    aes(x=1- log10(qglobal_cv.Normal), y=1-log10(qglobal_cv.Hypermutated), label=gene_name))+
    geom_point(data = dndsData[(dndsData$qglobal_cv.Normal >= plotThresh) & (dndsData$qglobal_cv.Hypermutated >= plotThresh),],
               aes(x=1- log10(qglobal_cv.Normal), y=1-log10(qglobal_cv.Hypermutated)))+
    
    geom_point(aes(x=1 - mean(log(dndsData[dndsData$qglobal_cv.Normal > 0,]$qglobal_cv.Normal)),
                   y=1 - mean(log(dndsData[dndsData$qglobal_cv.Hypermutated > 0,]$qglobal_cv.Hypermutated))), colour='orange', size=3)+
    #geom_text_repel(aes(x=1 - mean(log(dndsData[dndsData$qglobal_cv.Normal > 0,]$qglobal_cv.Normal)),
    #                    y=1 - mean(log(dndsData[dndsData$qglobal_cv.Hypermutated > 0,]$qglobal_cv.Hypermutated))), label='mean q_val', colour='orange')+
    
    xlim(minCoord,3)+
    ylim(minCoord,3)+
    geom_segment(aes(x=minCoord, xend=3, y=1- log10(.01), yend= 1- log10(.01)), colour='black', linetype=2)+
    geom_segment(aes(x=1- log10(.01), xend=1- log10(.01), y=minCoord, yend= 3), colour='black', linetype=2)+
    geom_segment(aes(x=minCoord, xend=3, y=1- log10(.1), yend= 1- log10(.1)), colour='black')+
    geom_segment(aes(x=1- log10(.1), xend=1- log10(.1), y=minCoord, yend= 3), colour='black')+
    ylab('1 minus q value in hypermutators')+
    xlab('1 minus q value in non-hypermutators')+
    ggtitle(title)
}

adjust_col_for_plot <- function(df){
  thresh <- 1e-10
  df$qHypermutated <- sapply(df$qglobal_cv.Hypermutated, function(x) if(1- log10(x) > 10){thresh}else{x})
  df$qNormal <- sapply(df$qglobal_cv.Normal, function(x) if(1- log10(x) > 10){thresh}else{x})
  return(df)
}

mergedEndo <- adjust_col_for_plot(mergedEndo)
mergedColo <- adjust_col_for_plot(mergedColo)
mergedGlioma <- adjust_col_for_plot(mergedGlioma)

fullDndsEndo <- make_dnds_plot(mergedEndo, 'Endometrial')
fullDndsColo <- make_dnds_plot(mergedColo, 'Colorectal')
fullDndsGlioma <- make_dnds_plot(mergedGlioma, 'Glioma')

weakDndsEndo <- make_dnds_weak_drivers_plot(mergedEndo, 'Endometrial')
weakDndsColo <- make_dnds_weak_drivers_plot(mergedColo, 'Colorectal')
weakDndsGlioma <- make_dnds_weak_drivers_plot(mergedGlioma, 'Glioma')

fullAligned = plot_grid(fullDndsEndo, fullDndsColo, fullDndsGlioma, ncol=3)
weakAligned = plot_grid(ggplot(), weakDndsEndo, weakDndsColo, weakDndsGlioma)

fullPlot = plot_grid(weakAligned, fullAligned, ncol=2, rel_widths=c(.5,1))

fullPlotWithTitle <- plot_grid(ggplot()+ggtitle('DNDS comparisson between Hypermutated and Normal')+theme(plot.title=element_text(hjust=.5)), 
                               fullPlot, ggplot() + labs(caption='runDnDsCv.R\nprepare_data_for_dnds_cv.ipynb'),
                               nrow=3, rel_heights = c(.1,1,.05))

ggsave('~/Desktop/plot.pdf', fullPlotWithTitle,  width = 14, height = 6, units = c("in"))


#
###
######
####
##

#QC 

write.table(dndsColorectalHyper$annotmuts, file='/Users/friedman/Desktop/hypermutationAnalysisProj/colorectalCancerHyperDNDS.tsv', quote=FALSE, sep='\t')
write.table(sel_cvEdnoH, file='/Users/friedman/Desktop/hypermutationAnalysisProj/endometrialCancerHyperDNDS.tsv', quote=FALSE, sep='\t')
write.table(sel_cvGliomaH, file='/Users/friedman/Desktop/hypermutationAnalysisProj/gliomaHyperDNDS.tsv', quote=FALSE, sep='\t')


#
####
########
##############
##################
########################
##################
##############
###########
#########
###

normalEndoGenes <- unique(mergedEndo[mergedEndo$qglobal_cv.Normal < .01, ]$gene_name)
hyperEndoGenes <- setdiff(unique(mergedEndo[mergedEndo$qglobal_cv.Hypermutated < .01, ]$gene_name), normalEndoGenes)
weakEndoGenes <- setdiff(unique(mergedEndo[(mergedEndo$qglobal_cv.Hypermutated < .1) | (mergedEndo$qglobal_cv.Normal < .1),]$gene_name)
                                , union(normalEndoGenes, hyperEndoGenes))

normalColoGenes <- unique(mergedColo[mergedColo$qglobal_cv.Normal < .01, ]$gene_name)
hyperColoGenes <- setdiff(unique(mergedColo[mergedColo$qglobal_cv.Hypermutated < .01, ]$gene_name), normalColoGenes)
weakColoGenes <- setdiff(unique(mergedColo[(mergedColo$qglobal_cv.Hypermutated < .1) | (mergedColo$qglobal_cv.Normal < .1),]$gene_name)
                         , union(normalColoGenes, hyperColoGenes))

normalGliomaGenes <- unique(mergedGlioma[mergedGlioma$qglobal_cv.Normal < .01, ]$gene_name)
hyperGliomaGenes <- setdiff(unique(mergedGlioma[mergedGlioma$qglobal_cv.Hypermutated < .01, ]$gene_name), normalGliomaGenes)
weakGliomaGenes <- setdiff(unique(mergedGlioma[(mergedGlioma$qglobal_cv.Hypermutated < .1) | (mergedGlioma$qglobal_cv.Normal < .1),]$gene_name)
                         , union(normalGliomaGenes, hyperGliomaGenes))




paste(weakGliomaGenes, collapse='","')

