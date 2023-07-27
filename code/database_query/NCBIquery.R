# Useful if i want to query NCBI DB 

rm(list = ls())

library(R.utils)

options(stringsAsFactors = F)

setwd("path to working dir")

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

###########################

filename_in <- paste0(dirDataset, "DEG.txt")
filename_out <- paste0(dirQuery, 'NCBI_gene_info_DEG.txt')

list <- read.table(".../code/Results/M0I/DEG.txt",
                   header = T, sep = "\t")
list <- data.frame(GeneSymbol = list$genes)

##################
fileURL <- "http://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"

destfile <- "Homo_sapiens_gene_info.txt.gz"

download.file(fileURL, destfile, method="auto")

gunzip(destfile)
##################

NCBI <- read.table("Homo_sapiens_gene_info.txt", sep ="\t", quote = "", header = T, comment.char = "", check.names = F)

df <- merge(list, NCBI, by.x = "GeneSymbol", by.y = "Symbol", all.x = T)
table(df$type_of_gene)
table(df$chromosome)

write.table(df, filename_out, sep ="\t", quote = F, row.names = FALSE)

file.remove("Homo_sapiens_gene_info.txt")
