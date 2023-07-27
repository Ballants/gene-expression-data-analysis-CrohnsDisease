## to install the package
#if(!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#
# BiocManager::install("GEOquery")
readr::local_edition(1)
library("GEOquery")

rm(list = ls())
setwd("path to working dir")

##### Set Results folder ####

dirRes <- "Results/"
if (!dir.exists(dirRes)){
  dir.create(dirRes)
}else{
  print(paste("The directory", dirRes, "already exists"))
}

dataset <- "M0I"

dirDataset <- paste0(dirRes,dataset, "/")
if (!dir.exists(dirDataset)){
  dir.create(dirDataset)
}else{
  print(paste("The directory", dirDataset, "already exists"))
}

###########################################
# STEP 1: Downloading data
###########################################

series = "GSE186582" # M0I
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186582
tmp <- getGEO(GEO = series) # it is a list
# show(tmp) 

set <- tmp$`GSE186582_series_matrix.txt.gz`
# phenotipic data
pData <- phenoData(set)

metadata_complete <- pData@data

metadata <- subset(metadata_complete, (`location:ch1` == "Ctrl" | `location:ch1` == "M0I") &
                  (`rutgeerts:ch1` == "Ctrl" | `rutgeerts:ch1` == "i3" | `rutgeerts:ch1` == "i4") &
                  (`rutgeertrec:ch1` == "Ctrl" | `rutgeertrec:ch1` == "Rec") &
                  (`postoperative anti tnf treatment:ch1` == "Ctrl" | `postoperative anti tnf treatment:ch1` == "No"))

aData <- assayData(set)
matrix <- data.frame(aData$exprs)
matrix<- 2^matrix

rm(pData, aData, tmp)

###########################################
# STEP 2: Preparing data
###########################################

annotation <- fData(set)
geneSymbol <- annotation$`Gene Symbol`
matrix <- matrix[annotation$ID,]
matrix <- aggregate(matrix, list(geneSymbol), "mean")
ind <- which(matrix$Group.1 == "")
if(length(ind) > 0) matrix <- matrix[-ind,]
rownames(matrix) <- matrix$Group.1
matrix <- matrix[, -1]
matrix[is.na(matrix)] <- 0

rm(set, annotation, ind, geneSymbol)

############################################
# STEP 3: Extracting case and control samples
###########################################

list <- split(metadata$geo_accession, metadata$`location:ch1`)
control <- list$Ctrl
case <- list$M0I

data <- matrix[,c(control, case)]
dataN <- matrix[,control]
dataC <- matrix[,case]

rm(list, matrix) 

############################################
#STEP 4: Export data
###########################################

write.table(data, paste0(dirDataset, "matrix.txt"),
            sep= "\t", col.names = NA, row.names = T, quote = F)
write.table(control, paste0(dirDataset, "control.txt"),
            sep= "\t", col.names = F, row.names = F, quote = F)
write.table(case, paste0(dirDataset, "case.txt"),
            sep= "\t", col.names = F, row.names = F, quote = F)
write.table(metadata_complete, paste0(dirDataset, "metadata_complete.txt"),
            sep= "\t", col.names = NA, row.names = T, quote = F)
write.table(metadata, paste0(dirDataset, "metadata.txt"),
            sep= "\t", col.names = T, row.names = F, quote = F)

# Remove downloaded data
unlink(series, recursive = TRUE)