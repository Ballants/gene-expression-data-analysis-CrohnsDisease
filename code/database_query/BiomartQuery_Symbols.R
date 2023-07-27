# source("https://bioconductor.org/biocLite.R")
# biocLite("biomaRt")
# https://bioconductor.org/packages/release/bioc/html/biomaRt.html
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("biomaRt")

rm(list=ls())

library(biomaRt)

setwd("path to working dir")
path_in <- ".../Results/"

options(stringsAsFactors=F)

################################################################

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
biomart_entrez_txt <- paste0(dirQuery, "biomart_entrez.txt")
biomart_geneBiotype_txt <- paste0(dirQuery, "biomart_geneBiotype.txt")
biomart_Ensembl_txt <- paste0(dirQuery, "biomart_Ensembl.txt") # it's a mess :)
biomart_refseq_txt <- paste0(dirQuery, "biomart_refseq.txt")

################################################################

# input
input_list <- read.table(".../Results/M0I/DEG.txt",
                         header = T, sep = "\t")
input_list <- data.frame(GeneSymbol = input_list$genes)


# download DB
ensembl <- useMart(biomart = "ensembl", dataset="hsapiens_gene_ensembl")
attributes = listAttributes(ensembl) # list of possibile output elements
filters = listFilters(ensembl) # list of possibile input elements

# download attributes

df_entrez <- getBM(attributes=c('hgnc_symbol','entrezgene_id'), 
                   filters='hgnc_symbol', values=input_list,
                   mart=ensembl)

df_Ensembl <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id', 'ensembl_transcript_id'),
                    filters='hgnc_symbol', values=input_list,
                    mart=ensembl)

# refseq id: 
# for each gene symbol there could be more than one refseq id, this happens because
# for each gene symbol there could be different isoforms that can give different 
# proteins
df_refseq <- getBM(attributes=c('hgnc_symbol','refseq_mrna'), 
                   filters='hgnc_symbol', values=input_list,
                   mart=ensembl)

ind <- which(df_refseq$refseq_mrna == "")
df_refseq <- df_refseq[-ind,]

df_gene_type <- getBM(attributes=c('hgnc_symbol','gene_biotype'), 
                      filters='hgnc_symbol', values=input_list,
                      mart=ensembl)
table(df_gene_type$gene_biotype) 

rm(attributes,filters,ensembl)

################################################################
# file output

write.table(df_entrez, file=biomart_entrez_txt, row.names=FALSE, 
            col.names=c("GeneSymbol","EntrezGene_ID"), sep = "\t", quote = FALSE)

write.table(df_gene_type, file=biomart_geneBiotype_txt, row.names=FALSE, 
            col.names=c("GeneSymbol","GeneBiotype"), sep = "\t", quote = FALSE)

write.table(df_Ensembl, file=biomart_Ensembl_txt, row.names=FALSE, 
            col.names=c("GeneSymbol","Ensembl_gene_ID", "Ensembl_transcript_ID"), sep = "\t", quote = FALSE)

write.table(df_refseq, file=biomart_refseq_txt, row.names=FALSE, 
            col.names=c("GeneSymbol","RefSeq"), sep = "\t", quote = FALSE)
