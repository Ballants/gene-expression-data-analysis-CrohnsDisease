################################################################
# Compute PCA of the differential expressed genes
################################################################

rm(list=ls())

#library(devtools)
#install_github("vqv/ggbiplot")

library(ggbiplot)

library(qcc)
library(ggpubr)
library(factoextra)
library(corrplot)
library(FactoMineR)
library(RColorBrewer)

setwd("path to working dir")
path_in <- ".../Results/"

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

# output files
file_score_plot <- paste0(dirDataset, "score_plot.pdf")
file_pareto_scree_plot <- paste0(dirDataset, "pareto_scree_plot.pdf")
file_loading_plot <- paste0(dirDataset, "loading_plot.pdf") # it's a mess :)
file_contribution_plot <- paste0(dirDataset, "PC_contribution_plot.pdf")

################################################################
# 1. Importing data

data <- read.table(paste0(path_in, dataset, "/matrix_DEG.txt"),
                   header = T, sep = "\t", quote = "",
                   check.names = F, row.names = 1)

list_control <- read.table(paste0(path_in, dataset, "/control.txt"),
                          header = F, sep = "\t", quote = "", check.names = F)$V1
list_case <- read.table(paste0(path_in, dataset, "/case.txt"),
                         header = F, sep = "\t", quote = "", check.names = F)$V1

data <- t(data[,c(list_case, list_control)])

groups <- c(rep("case", length(list_case)), rep("control", length(list_control)))

################################################################
# 2. Apply PCA
# Rows of data correspond to observations (samples), columns to variable (genes)

pca <- prcomp(data, center = T, scale = T, retx = T)

################################################################
# 3. Compute score and store score plot
# (scores = the coordinates of old data (observations) in the new)

# pca$x = t(data)*pca$rotation
scores <- pca$x

pdf(file_score_plot, width = 5, height = 5)
g <- ggbiplot(pca, obs.scale = 1, var.axes = F,
              ellipse = T, groups = groups)
print(g)
dev.off()

################################################################
# 4. Compute eigenvalue
# eigenvalues of the covariance matrix ordered in decreasing order
# (from the largest to the smaller)
eigenvalue = pca$sdev^2

# variance explained by each PC
varS <- round(eigenvalue/sum(eigenvalue)*100, 2)
names(varS) = paste0('PC', seq(1, length(varS)))

pdf(file_pareto_scree_plot, width = 5, height = 5)

# pareto chart
# in which u have all the cumulative variance of all the PC
pareto.chart(varS[1:10])

# scree plot
fviz_eig(pca, addlabels = TRUE)
# mean_lambda = 0.7*mean(eigenvalue) # threshold for selecting PCs
dev.off()

################################################################
# 5. Compute loadings (coefficients of each PC)
# the matrix of variable loadings (a matrix whose columns contain the eigenvectors)

loadings <- pca$rotation

# Contribution of variables to the PCs
contrib_var <- get_pca_var(pca)$contrib
colnames(contrib_var) <- paste0('PC', seq(1, ncol(contrib_var)))

contrib_var <- contrib_var[order(contrib_var[,"PC1"], decreasing = T),]

# can be skipped? 
corrplot(contrib_var[1:5, 1:5], is.corr = FALSE,
         tl.col = "black", method = "circle",
         col = brewer.pal(n = 9, name = "BuPu"),
         addCoef.col = "black")

# this must be done
pdf(file_contribution_plot, width = 5, height = 5)

# to PC1
fviz_contrib(pca, choice = "var", axes = 1, top = 10)
# to PC2
fviz_contrib(pca, choice = "var", axes = 2, top = 10)
# to PC3
fviz_contrib(pca, choice = "var", axes = 3, top = 10)

# total (PC1 + PC2 + PC3)
fviz_contrib(pca, choice = "var", axes = 1:3, top = 10)

dev.off()

################################################################

