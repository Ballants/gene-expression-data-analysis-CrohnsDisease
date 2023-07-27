################################################################
# Functional enrichment analysis - enrichR
################################################################

#install.packages("devtools") 
#devtools::install_github("wjawaid/enrichR")

rm(list=ls())

library(enrichR)
library(ggplot2)
library(forcats)
library(stringr)

################################################
setwd("path to working dir")
source(".../code/database_query/enrichR_UP_DOWN/getEnrichment.R")
source(".../code/database_query/enrichR_UP_DOWN/getEnrichmentPlot.R")
################################################

dataset <- "M0I"

dirRes <- "Results/"

if (!dir.exists(dirRes)){
  dir.create(dirRes)
}else{
  print(paste("The directory", dirRes, "already exists"))
}

dirDataset <- paste0(dirRes, dataset, "/")

if (!dir.exists(dirDataset)){
  dir.create(dirDataset)
}else{
  print(paste("The directory", dirDataset, "already exists"))
}

dirEnrich <- paste0(dirDataset, "Functional_Enrichment/")

if (!dir.exists(dirEnrich)){
  dir.create(dirEnrich)
}else{
  print(paste("The directory", dirEnrich, "already exists"))
}

################################################
top_term <- 10
thr_pval <- 0.05
################################################
file_input_list <- paste0(dirDataset, "DEG.txt")

# websiteLive <- getOption("enrichR.live")
# if (websiteLive) dbs <- listEnrichrDbs() #list of all dbs

# TODO modify db

# Transcription tab: 
#   to see if the genes can be regulated by some "external" transcription factors,
#   if the protein region of the genes could bind some transcription factor.
# TRANSFAC_and_JASPAR_PWMs

# Pathways tab:
# Reactome_2022
# KEGG_2021_Human

# Ontologies tab:
# GO_Cellular_Component_2021: where the genes are located in the cell
# GO_Molecular_Function_2021: which are the functions that our the genes can perform 
# GO_Biological_Process_2021: in which processes the genes are involved

# Disease/Drugs tab:
# DisGeNET
# GWAS_Catalog_2023

dbs <- c("TRANSFAC_and_JASPAR_PWMs", 
         "KEGG_2021_Human", "Reactome_2022",
         "GO_Cellular_Component_2023", "GO_Molecular_Function_2023", "GO_Biological_Process_2023", 
         "DisGeNET", "GWAS_Catalog_2023")

input_list <- read.table(file_input_list, sep = "\t", header = T, check.names = F, quote = "")
list <- split(input_list$genes, input_list$direction)

df <- lapply(list, function(x){
  enrichr(x, dbs)
})



getEnrichment(df$UP, "UP")
getEnrichment(df$DOWN, "DOWN")



################################################



# Enrichment analysis of miRNA targets

genes <- c("CYP1B1", "CTGF", "LAMC1", "LITAF", "NID1", "ELK3", "CCL2", "PEA15", "TUBB6", "ZCCHC24", "NFKBIZ", "VNN2", "LAMA4",
                 "PRRX1", "COL4A1", "IL6", "SERPINE1", "RBMS1", "MTCL1", "MSRB3", "FAM129A", "CD55", "GNB4", "ENAH", "DNM3", "EHD2",
                 "TMEM45A", "MCAM", "RBPMS", "GBP1", "PLA2G7", "CKAP4", "ADAMTS1", "BCL6", "CEMIP", "PTGS2", "COL1A1", "KANK2", 
                 "MUC1", "FAM83D", "NR4A3", "EGR2", "CTHRC1", "GREM1", "C4BPB", "RIPOR3", "RAB31", "FSTL1", "TUBA1A", "CLMP", "NCF2",
                 "SAMSN1", "RNF144B", "CDH11", "TRPS1", "AHNAK2", "PTGS1", "KLF2", "CAVIN3", "CAVIN1", "CXCL8", "CXCL1", "CYR61", 
                 "IRAK3", "MMP19", "EGR1", "AMOTL1", "CD274", "PTPN14")

dbs <- c("KEGG_2021_Human", "Reactome_2022",
         "DisGeNET")

df <- enrichr(genes, dbs)

KEGG <- df$KEGG_2021_Human
Reactome <- df$Reactome_2022
DisGeNET <- df$DisGeNET

getEnrichment(df, "miRNA")