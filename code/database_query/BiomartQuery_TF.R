# In EnrichR u can see the transcription activity, but what u see in the "transcription" factor
# section is if your list of genes can be regulated by some "external" transcription factors,
# if the protein region of your genes could bind some transcription factor.

# IN THIS SCRIPT, WE ARE LOOKING FOR GENES (in our list of dysregulated genes) THAT COULD BE TRANSCRIPTION FACTORS.

rm(list = ls())

library(biomaRt)

setwd("path to working dir")

options(stringsAsFactors=F)

###########################

dirRes <- "Results/"

if (!dir.exists(dirRes)){
  dir.create(dirRes)
}else{
  print(paste("The directory", dirRes, "already exists"))
}

dataset <- "M0I"

dirDataset <- paste0(dirRes, dataset, "/")
if (!dir.exists(dirDataset)){
  dir.create(dirDataset)
}else{
  print(paste("The directory", dirDataset, "already exists"))
}

dirQuery <- paste0(dirDataset, "Other database query/")

if (!dir.exists(dirQuery)){
  dir.create(dirQuery)
}else{
  print(paste("The directory", dirQuery, "already exists"))
}

# output files
TF_file <- paste0(dirQuery, "TF.txt")

###########################

list <- read.table(".../code/Results/M0I/DEG.txt",
                   header = T, sep = "\t")
list <- data.frame(GeneSymbol = list$genes)

# download database
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
attributes = listAttributes(ensembl) # list of possible output elements
filters = listFilters(ensembl) # list of possible input elements

# download attributes
df <- getBM(attributes = c('hgnc_symbol', 'name_1006', "namespace_1003", 'go_id'), 
            filters = 'hgnc_symbol', values = list, mart = ensembl)

colnames(df) <- c("GeneSymbol", "GO_term", "GO domain", "GO ID")

rm(attributes,filters,ensembl)

###########################

# search for TF
ind_TF <- grep("transcription factor activity", df$GO_term)

TF <- unique(df[ind_TF, "GeneSymbol"])

###########################

# write file output
write.table(TF, file=TF_file, col.names = F, row.names = FALSE, quote = FALSE)
