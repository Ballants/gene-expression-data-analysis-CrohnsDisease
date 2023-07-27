################################################################
# Gene expression data analysis
################################################################# 
# Aim: finding a list of genes that are differential expressed between case
# and control samples and characterize them.

rm(list = ls())
# when you import character in R they are seen as factors
options(stringsAsFactors = F)  
setwd("path to working dir")

################################################################

library(stringr) # used to import the function "str_extract"
library(pheatmap) # used to plot the heatmap

################################################################

dirRes <- "Results/"

if (!dir.exists(dirRes)){
  dir.create(dirRes)
}else{
  print(paste("The directory", dirRes, "already exists"))
}

# the name of the case that we are studying
dataset <- "M0I"

dirDataset <- paste0(dirRes, dataset, "/")
if (!dir.exists(dirDataset)){
  dir.create(dirDataset)
}else{
  print(paste("The directory", dirDataset, "already exists"))
}

################################################################

# parameters setting - depend from the chosen dataset
prc_IQR <- 0.2
thr_fc <- 2
thr_pval <- 0.05
paired <- FALSE

# output 
filename_DEG <- paste0(dirDataset, "DEG.txt")
filename_matrix_DEG <- paste0(dirDataset, "matrix_DEG.txt")
filename_heatmap <- paste0(dirDataset, "heatmap.pdf")

################################################################
# STEP 1: Importing data
################################################################

data <- read.table(paste0(dirDataset, "matrix.txt"), header = T, sep = "\t",
                   quote = "", check.names = F, row.names = 1)

list_control <- read.table(paste0(dirDataset, "control.txt"), header = F,
                          quote = "", check.names = F)$V1
list_case <- read.table(paste0(dirDataset, "case.txt"), header = F,
                         sep = "\t", quote = "", check.names = F)$V1

dataN <- data[, list_control]
dataC <- data[, list_case]
data <- data[, c(list_control, list_case)]

metadata <- read.table(paste0(dirDataset, "metadata.txt"), header = T, 
                       sep = "\t", quote = "", check.names = F)

################################################################
# STEP 2: Analysis
################################################################

genes <- row.names(data)
pz <- colnames(data)
pzN <- colnames(dataN)
pzC <- colnames(dataC)

#############################
# Pre-processing:
#############################

##############
# STEP 2.0: For each matrix (data_control and data_case), average the
#           expression data for each gene and remove those ones whose 
#           mean is equal to zero in the two matrices simultaneously.

# General step: compute a variable (ex: mean, iqr, fc, p_value),
#               choose a threshold value,
#               filter data for the given threshold,
#               remove results from all the matrices we have.

overall_mean <- rowMeans(data)
ind <- which(overall_mean == 0)
if(length(ind) > 0){
  dataN <- dataN[-ind,]
  dataC <- dataC[-ind,]
  data <- data[-ind,]
  genes <- genes[-ind]
}

rm(ind)

##############
# STEP 2.1: Logarithmic transformation. 
#           Make the log 2 of the (data + 1)
dataN <- log2(dataN + 1)
dataC <- log2(dataC + 1)
data <- log2(data + 1)

##############
# STEP 2.2: calculate IQR value for each gene and the 10th percentile of
#           the IQR distribution (choose an appropriate threshold to
#           eliminate those ones that slightly vary: for instance, remove
#           those genes with IQR <10th (or higher) percentile according
#           to the threshold suggested by the graph of the frequency of
#           IQR distribution)

variation <- apply(data, 1, IQR)

# IQR filtering
thr_pcr <- quantile(variation, prc_IQR)
thr_pcr

#par(mar = c(2, 2, 2, 2))  # Adjust the margin sizes to smaller values

# In order to see how many genes we are filtering out we can perform an histogram:
hist(variation,
     main = "IQR frequency distribution", breaks = 100,
     xlab = "IQR value", ylab = "Frequency", col = "blue")
abline(v = thr_pcr, lty = 2, lwd = 4, col = "red")


ind <- which(variation <= thr_pcr)
if(length(ind) > 0){
  dataN <- dataN[-ind,]
  dataC <- dataC[-ind,]
  data <- data[-ind,]
  genes <- genes[-ind]
}

rm(ind)

#############################
# Filtering:
#############################

##############
# STEP 2.3: Calculate the fold-change (case/control, data are on log scale!
# => logFC = mean(case) - mean(control))

logFC <- rowMeans(dataC) - rowMeans(dataN)

hist(logFC,
     main = "FC (logarithmic) frequency distribution", breaks = 100,
     xlab = "log FC", ylab = "Frequency", col = "blue")
abline(v = c(-log2(thr_fc), log2(thr_fc)), lty = 2, lwd = 4, col = "red")

##############
# STEP 2.4: Set a fold-change threshold and remove all transcripts 
# with |logFC|<log(threshold).
# Filter out the middle part of our histogram:

ind <- which(abs(logFC) < log2(thr_fc))

if(length(ind)>0){
  dataN <- dataN[-ind,]
  dataC <- dataC[-ind,]
  data <- data[-ind,]
  genes <- genes[-ind]
  logFC <- logFC[-ind]
}

rm(ind)

##############
# STEP 2.5: If the samples size is n>30, you can exploit the hypothesis that the two distributions 
#   (cancer tissues and control tissues) are control and perform a paired Student’s t-test:
#     − Null Hypothesis: the two samples are from the same distribution (the averages are the same)
#     − Alternative Hypothesis: the two samples come from two different distributions (the averages are different)

N <- ncol(dataN)
M <-  ncol(dataC)
pval <- apply(data, 1, function(x){
  t.test(x[1:N], x[(N+1):(M+N)], paired=paired)$p.value
})

##############
# STEP 2.6: Adjust the p-values (e.g., by using the False Discovery Rate)
pval_adj <- p.adjust(pval, method='fdr')

# Draw a vulcano plot
plot(logFC, -log10(pval_adj),
     main = "Volcano plot",
     xlim = c(-8, 8),
     ylim = c(0, 60),
     xlab = "log2 fold change",
     ylab = "-log10 p-value")
abline(v = c(-log2(thr_fc), log2(thr_fc)), lty = 2, lwd = 4, col = "blue")
abline(h = -log10(thr_pval), lty = 2, lwd = 4, col = "red")


# Remove those genes with adjusted-pvalue greater than the
# chosen threshold of significance level (e.g., 0.05 or 0.01)
ind <- which(pval_adj > thr_pval)
if(length(ind)>0){
  dataN <- dataN[-ind,]
  dataC <- dataC[-ind,]
  data <- data[-ind,]
  genes <- genes[-ind]
  logFC <- logFC[-ind]
  pval <- pval[-ind]
  pval_adj <- pval_adj[-ind]
}

rm(ind)

################################################################
# STEP 3: Exporting results
################################################################
# Save the obtained results in a txt file reporting, for each differential
# expressed gene, the gene name, the adjusted p-value, and the logFC.

# adding a column to specify if the genes is up or down regulated
direction <- ifelse(logFC>0, "UP", "DOWN")

results <- data.frame(genes = genes,
                      pvalue = pval, pval_adj = pval_adj,
                      logFC = logFC, direction = direction)
# ordering the results
results <- results[order(results$logFC, decreasing = T),]

write.table(results, file = filename_DEG, row.names = F,
            sep = "\t", quote = F)

write.table(data, file = filename_matrix_DEG, row.names = T, 
            col.names = NA, sep = "\t", quote = F)

#write.table(pzN_com, file = filename_list_control, row.names = F, 
#            col.names = F, sep = "\t", quote = F)

#write.table(pzC_com, file = filename_list_case, row.names = F, 
#            col.names = F, sep = "\t", quote = F)

################################################################
# STEP 4: PLOTs
################################################################

####################
# Volcano plot

# Draw another vulcano plot with remaining data after the filtering

plot(logFC, -log10(pval_adj),
     main = "Volcano plot",
     xlim = c(-8, 8),
     ylim = c(0, 60),
     xlab = "log2 fold change",
     ylab = "-log10 p-value")

abline(h = -log10(thr_pval), lty = 2, lwd = 4, col = "red")
abline(v = c(-log2(thr_fc), log2(thr_fc)), lty = 2, lwd = 4, col = "blue")

####################
# Box plot

# Draw boxplot of a selected differential expressed gene between cancer and 
# control tissues and report on the plot its corresponding p-value

ind <- which.max(results$logFC) # most up-regulated gene
gene_id <- genes[ind]

boxplot(as.numeric(dataN[ind,]),
        as.numeric(dataC[ind,]),
        main = paste0(gene_id, ",", " adjusted p_value = ", 
                      format(pval_adj[ind], digits = 2)),
        notch = T,
        ylab = "Gene expression value",
        xlab = "Condition",
        col = c("green", "orange"),
        names = c("control", "case"),
        pars = list(boxwex = 0.3, staplewex = 0.6))

ind <- which.min(results$logFC) # most down-regulated gene
gene_id <- genes[ind]

boxplot(as.numeric(dataN[ind,]),
        as.numeric(dataC[ind,]),
        main = paste0(gene_id, ",", " adjusted p_value = ", 
                      format(pval_adj[ind], digits = 2)),
        notch = T,
        ylab = "Gene expression value",
        xlab = "Condition",
        col = c("green", "orange"),
        names = c("control", "case"),
        pars = list(boxwex = 0.3, staplewex = 0.6))


ind <- 200 # random gene
gene_id <- genes[ind]

boxplot(as.numeric(dataN[ind,]),
        as.numeric(dataC[ind,]),
        main = paste0(gene_id, ",", " adjusted p_value = ", 
                      format(pval_adj[ind], digits = 2)),
        notch = T,
        ylab = "Gene expression value",
        xlab = "Condition",
        col = c("green", "orange"),
        names = c("control", "case"),
        pars = list(boxwex = 0.3, staplewex = 0.6))

####################
# Pie chart

count <- table(results$direction)
pie(count,
    label = paste0(names(count), " ",
                   round(100 * count/sum(count), 2), "%"),
    col = c("blue", "gold"))    

####################
# Clusterize the differential expressed genes and samples and
# draw the heatmap of the differential expressed genes

# create a file of annotations
annotation <- data.frame(condition = metadata$`location:ch1`)

# row names: id (those which compare in your matrix)
rownames(annotation) <- metadata$geo_accession 

vect_color <- c("green", "violet")
names(vect_color) <- unique(annotation$condition) 

annotation_color <- list(condition = vect_color)

heatmap_result <- pheatmap(data, scale = "row",
         border_color = NA,
         cluster_cols = T, cluster_rows = T,
         clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation", 
         clustering_method = "average", 
         annotation_col = annotation,
         annotation_colors = annotation_color,
         color = colorRampPalette(colors = 
                                    c("blue", "blue3", "black",
                                            "yellow3", "yellow"))(100),
                                            show_rownames = F, show_colnames = F,
         cutree_cols = 2, cutree_rows = 2,
         width = 10, height = 10,
         filename = filename_heatmap,
         keep.dendro = T
)