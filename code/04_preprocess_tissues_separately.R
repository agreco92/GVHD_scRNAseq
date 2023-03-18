# libraries ---------------------------------------------------------------
set.seed(210201)

source('./bin/useful_functions.R')
library(scater)
library(scran)
library(Seurat)
library(org.Mm.eg.db)

library(tidyverse)
library(ggthemes)
library(ggpubr)
library(gridExtra)

# import global objects ------------------------------------------------------------------
output_folder = './data/1_preprocessing/'
dir.create(output_folder, recursive = T)

sce <- readRDS('./data/1_preprocessing/02_all_tissues_hires.rds')
sct <- readRDS('./data/1_preprocessing/01_sctransform_residuals.rds')
hvgs <- readRDS('./data/1_preprocessing/01_sctransform_hvgs.rds')

tissues = unique(sce$tissue_hires)
# exclude LP_IEL as redundant with sorted IELs
tissues = tissues[tissues != 'LP_IEL']

# tissue-wise dimensionality reduction and clusters  ------------------------------------------
# 
for(ti in tissues){
  tissue_sce <- sce[,sce$tissue_hires == ti]
  tissue_sce_sct <- SingleCellExperiment(
    assays = list(counts = sct[hvgs,sce$tissue_hires == ti], 
                  logcounts = sct[hvgs,sce$tissue_hires == ti]))
  tissue_sce_sct <- runPCA(tissue_sce_sct, subset_row = hvgs)
  tissue_sce_sct <- runUMAP(tissue_sce_sct, dimred = 'PCA', n_dimred = 15, min_dist = 0.5)
  reducedDims(tissue_sce) <- reducedDims(tissue_sce_sct)
  
  # clustering
  for(res in c(0.1,0.2,0.25,0.3, 0.4, 0.5, 0.6, 0.75)){
    tissue_sce <- seurat_clusters(tissue_sce, dim.red = 'PCA', ndimred = 30, resolution_param = res, 
                                  colname = paste0('louvain', res))}
  
  # remove clustering solutions containing only one cluster
  cluster_cols <- grep('^louvain', colnames(colData(tissue_sce)), value = T)
  degenerate_clusters <- cluster_cols[which(sapply(cluster_cols, function(x){dim(table(colData(tissue_sce)[[x]]))}) < 2)]
  colData(tissue_sce)[,degenerate_clusters] <- NULL
  cluster_cols <- cluster_cols[!cluster_cols %in% degenerate_clusters]
  
  # save sce objects
  saveRDS(tissue_sce, 
          paste0('./data/2_tissue_preprocessing/',ti, '.rds'))
}

