# includes ----------------------------------------------------------------
set.seed(210201)
library(scater)
library(scran)
library(tidyverse)
library(ggthemes)
library(ggrepel)
library(ggrastr)
library(patchwork)
source('./bin/useful_functions.R')

sce = readRDS('./data/1_preprocessing/02_all_tissues_hires.rds')
dissoc_genes <- read.csv('./data/gene_sets/dissociation_genes_full.txt', skip = 1)[,1]

# compute tissue markers pairwise -------------------------------------------------------------------
tissue_names <- unique(sce$tissue)

combos <- combn(tissue_names,2) # pairwise combinations of tissues
n_markers = 10

# create list with tissues 
tissue_list = list()
for(t in tissue_names){
  tissue_list[[t]] = sce[,sce$tissue ==t]
}

top_hits <- c()
n_markers = 10
for (i in 1:ncol(combos)){
  t1 <- tissue_list[[combos[1,i]]]
  t2 <- tissue_list[[combos[2,i]]]
  
  # remove lowly expressed genes
  t1_pos <- rowSums(sign(counts(t1)))/ncol(t1)
  t2_pos <- rowSums(sign(counts(t2)))/ncol(t2)
  low_exp <- t1_pos < 0.25 & t2_pos <0.25 
  dissoc_genes_index <- rownames(t1) %in% dissoc_genes
  
  # create SCE with tissue pair
  t <- SingleCellExperiment(
    assays = list(counts = cbind(counts(t1[!low_exp & !dissoc_genes_index,]), 
                                 counts(t2[!low_exp & !dissoc_genes_index]))))
  t$tissue = c(rep(combos[1,i], ncol(t1)), rep(combos[2,i], ncol(t2)))
  
  # compute tissue markers between t1, t2
  markers <- findMarkers(t, groups = t$tissue, test.type = 'binom',
                         assay.type = 'counts', direction = 'up')
  
  # store top 10 markers
  top_hits <- c(top_hits, rownames(markers[[1]][1:n_markers,]), 
                rownames(markers[[2]][1:n_markers,]))
}

# plot tissue markers with detected at least twice
top_hits <- table(top_hits)
tissue_markers <- names(top_hits[top_hits>2])

# create dataframe with % of expressing cells for each markers
tm_df <- data.frame(row.names = tissue_markers)
tm_df[,tissue_names] <- NA

for(t in tissue_names){
  obj <- tissue_list[[t]]
  for(m in tissue_markers){
    # positive cells
    tm_df[m,t] <- rowSums(sign(counts(obj[m,])))/ncol(obj)
  }
}
cols <- viridis::plasma(101)
pheatmap::pheatmap(tm_df, 
                   filename = './code/figures/6b_tissue_markers.pdf', 
                   color = cols, border_color = NA)

# compute tissue markers pairwise (split colon) -------------------------------------------------------------------
tissue_names <- unique(sce$tissue_hires)
tissue_names <- tissue_names[tissue_names!= 'LP_IEL']
combos <- combn(tissue_names,2) # pairwise combinations of tissues
n_markers = 7

# create list with tissues 
tissue_list = list()
for(t in tissue_names){
  tissue_list[[t]] = sce[,sce$tissue_hires ==t]
}

top_hits <- c()
n_markers = 7
for (i in 1:ncol(combos)){
  t1 <- tissue_list[[combos[1,i]]]
  t2 <- tissue_list[[combos[2,i]]]
  
  # remove lowly expressed genes
  t1_pos <- rowSums(sign(counts(t1)))/ncol(t1)
  t2_pos <- rowSums(sign(counts(t2)))/ncol(t2)
  low_exp <- t1_pos < 0.25 & t2_pos <0.25 
  dissoc_genes_index <- rownames(t1) %in% dissoc_genes
  
  # create SCE with tissue pair
  t <- SingleCellExperiment(
    assays = list(counts = cbind(counts(t1[!low_exp & !dissoc_genes_index,]), 
                                 counts(t2[!low_exp & !dissoc_genes_index]))))
  t$tissue = c(rep(combos[1,i], ncol(t1)), rep(combos[2,i], ncol(t2)))
  
  # compute tissue markers between t1, t2
  markers <- findMarkers(t, groups = t$tissue, test.type = 'binom',
                         assay.type = 'counts', direction = 'up')
  
  # store top 10 markers
  top_hits <- c(top_hits, rownames(markers[[1]][1:n_markers,]), 
                rownames(markers[[2]][1:n_markers,]))
}

# plot tissue markers with detected at least twice
top_hits <- table(top_hits)
tissue_markers <- names(top_hits[top_hits>2])

# create dataframe with % of expressing cells for each markers
tm_df <- data.frame(row.names = tissue_markers)
tm_df[,tissue_names] <- NA

for(t in tissue_names){
  obj <- tissue_list[[t]]
  for(m in tissue_markers){
    # positive cells
    tm_df[m,t] <- rowSums(sign(counts(obj[m,])))/ncol(obj)
  }
}

cols <- viridis::plasma(101)
pheatmap::pheatmap(tm_df, 
                   filename = './code/figures/s6g_tissue_markers_colonsplit.pdf', 
                   color = cols, border_color = NA)

