# include -------------------------------------------------------------------------

set.seed(200907)
# 
library(scater)
library(scran)
library(Seurat)
# 
library(tidyverse)
library(ggthemes)
library(ggpubr)
source('./bin/useful_functions.R')

# import skin -------------------------------------------------------------

skin <- readRDS('./data/2_tissue_preprocessing/Skin.rds')
dissociation_genes <- as.character(read.csv('./data/dissociation_genes.txt', header =F)[,1])
dissociation_genes <- intersect(dissociation_genes, rownames(skin))

# score dissociation signature
seu <- as.Seurat(skin)
seu <- AddModuleScore(seu, features = list(dissociation_genes), name = 'dissociation_score')
skin$dissociation_score <- seu@meta.data$dissociation_score1
plotColData(skin, y = 'dissociation_score', x = 'louvain0.5')

# filter out dissociation cluster
skin <- skin[,skin$louvain0.5!= 3]


# repeat dr,  clustering ----------------------------------------------------------------------
sct <- readRDS('./data/1_preprocessing/01_sctransform_residuals.rds')
hvgs <- readRDS('./data/1_preprocessing/01_sctransform_hvgs.rds')

skin_sct <- sct[,colnames(skin)]

skin_sct <- SingleCellExperiment(list(counts = skin_sct[hvgs,], 
                                            logcounts = skin_sct[hvgs,]))
skin_sct <- runPCA(skin_sct, subset_row = hvgs)
skin_sct <- runUMAP(skin_sct, dimred = 'PCA',n_dimred = 15, min_dist = 0.5)
reducedDims(skin) <- reducedDims(skin_sct)


# repeat clustering
for(res in c(0.1,0.2,0.25,0.3, 0.4, 0.5, 0.6, 0.75)){
  skin <- seurat_clusters(skin, dim.red = 'PCA', ndimred = 30, resolution_param = res, 
                                colname = paste0('louvain', res))}

# remove clustering solutions containing only one cluster
cluster_cols <- grep('^louvain', colnames(colData(skin)), value = T)
degenerate_clusters <- cluster_cols[which(sapply(cluster_cols, function(x){dim(table(colData(skin)[[x]]))}) < 2)]
colData(skin)[,degenerate_clusters] <- NULL
cluster_cols <- cluster_cols[!cluster_cols %in% degenerate_clusters]

# save ----------------------------------------------------------------------------------------

# save skin object
saveRDS(skin, './data/2_tissue_preprocessing/Skin_filtered.rds')
