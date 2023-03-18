set.seed(210201)

library(scater)
library(scran)
library(ggthemes)
library(ggrepel)
library(patchwork)
library(tidyverse)
library(org.Mm.eg.db)
source('./bin/useful_functions.R')

spleen <- readRDS('./data/2_tissue_preprocessing/Spleen.rds')
skin <- readRDS('./data/2_tissue_preprocessing/Skin_filtered.rds')
lp <- readRDS('./data/2_tissue_preprocessing/LP.rds')
colon_lp <- readRDS('./data/2_tissue_preprocessing/Colon_LP.rds')
colon_iel <- readRDS('./data/2_tissue_preprocessing/Colon_IEL.rds')
iel <- readRDS('./data/2_tissue_preprocessing/IEL.rds')


# hvg detection -----------------------------------------------------------
# marker detection using sct is strongly biased towards lowly expressed markers.
# For this reason, marker overlap is computed on on the 
# union of all scran-based hvgs
tissue_names = c('spleen', 'colon_lp', 'colon_iel', 'iel', 'lp', 'skin')

all_hvgs = lapply(X = tissue_names,
       function(t){
         tissue_sce <- get(t)
         genevar <- modelGeneVar(tissue_sce)
         hvg <- getTopHVGs(genevar, n = 1000)
         return(hvg)
       }) %>% Reduce(f = dplyr::union)

available_genes = 
  Reduce(
    f = intersect,
    x = lapply(X = tissue_names, 
               FUN = function(t){
                 tissue_sce <- get(t)
                 return(rownames(tissue_sce))}))

selected_features <- intersect(available_genes, all_hvgs)


# marker_detection --------------------------------------------------------

for(t in tissue_names){
  assign(t, get(t)[selected_features,])}

#### compute markers #####
# markers
max_markers <- 100
min_prop <- 0.5
fdr_threshold <- 0.01 
spleen_res <- 'louvain0.5'
lp_res <- 'louvain0.4'
iel_res <- 'louvain0.4'
coloniel_res <- 'louvain0.75'
colonlp_res <- 'louvain0.4'
skin_res <- 'louvain0.3'

# markers identification --------------------------------------------------
# colon_lp
colon_lp_markers <- findMarkers(colon_lp, groups = colData(colon_lp)[,colonlp_res], direction = 'up', 
                                test.type = 'wilcox', pval.type = 'some', min.prop = min_prop, 
                                subset.row = selected_features)

# store markers for each cluster in a matrix. When the number of markers is lower than 100,
# the column is filled with "dummy" values that will be discarded later
colon_lp_markers <- sapply(colon_lp_markers, function(x){
  tmp <- rownames(x)
  tmp[x$FDR >= fdr_threshold] <- 'dummy'
  tmp <- tmp[1:max_markers]
  return(tmp)
})
colnames(colon_lp_markers) <- paste0('colon_lp_', colnames(colon_lp_markers))

# colon_iel
colon_iel_markers <- findMarkers(colon_iel, groups = colData(colon_iel)[,coloniel_res], direction = 'up', 
                                 test.type = 'wilcox', pval.type = 'some', min.prop = min_prop, 
                                 subset.row = selected_features)
# low cell numbers require looser threshold for marker detection,
# thus the pval.type parameter is set to "all"
colon_iel_markers <- findMarkers(colon_iel, groups = colData(colon_iel)[,coloniel_res], direction = 'up', 
                                 test.type = 'wilcox',  
                                 subset.row = selected_features)

colon_iel_markers <- sapply(colon_iel_markers, function(x){
  tmp <- rownames(x)
  tmp[x$FDR >= fdr_threshold] <- 'dummy'
  tmp <- tmp[1:max_markers]
  return(tmp)
})
colnames(colon_iel_markers) <- paste0('colon_iel_', colnames(colon_iel_markers))

# lp
lp_markers <- findMarkers(lp, groups = colData(lp)[,lp_res], direction = 'up', 
                          test.type = 'wilcox', pval.type = 'some', min.prop = min_prop, 
                          subset.row = selected_features)

lp_markers <- sapply(lp_markers, function(x){
  tmp <- rownames(x)
  tmp[x$FDR >= fdr_threshold] <- 'dummy'
  tmp <- tmp[1:max_markers]
  return(tmp)
})
colnames(lp_markers) <- paste0('lp_', colnames(lp_markers))


# iel
iel_markers <- findMarkers(iel, groups = colData(iel)[,iel_res], direction = 'up', 
                           test.type = 'wilcox', pval.type = 'some', min.prop = min_prop, 
                           subset.row = selected_features)

iel_markers <- sapply(iel_markers, function(x){
  tmp <- rownames(x)
  tmp[x$FDR >= fdr_threshold] <- 'dummy'
  tmp <- tmp[1:max_markers]
  return(tmp)
})
colnames(iel_markers) <- paste0('iel_', colnames(iel_markers))

# skin
skin_markers <- findMarkers(skin, groups = colData(skin)[,skin_res], direction = 'up', 
                            test.type = 'wilcox', pval.type = 'some', min.prop = min_prop, 
                            subset.row = selected_features)

skin_markers <- sapply(skin_markers, function(x){
  tmp <- rownames(x)
  tmp[x$FDR >= fdr_threshold] <- 'dummy'
  tmp <- tmp[1:max_markers]
  return(tmp)
})
colnames(skin_markers) <- paste0('skin_', colnames(skin_markers))

# spleen
spleen_markers <- findMarkers(spleen, groups = colData(spleen)[,spleen_res], direction = 'up', 
                              test.type = 'wilcox', pval.type = 'some', min.prop = min_prop, 
                              subset.row = selected_features)

spleen_markers <- sapply(spleen_markers, function(x){
  tmp <- rownames(x)
  tmp[x$FDR >= fdr_threshold] <- 'dummy'
  tmp <- tmp[1:max_markers]
  return(tmp)
})
colnames(spleen_markers) <- paste0('spleen_', colnames(spleen_markers))

# compute overlap between markers -------------------------------------------------------

# create a matrix with features as rows and clusters as columns -------------------------------

all_markers <- Reduce(cbind, list(spleen_markers, colon_lp_markers, 
                                   colon_iel_markers,
                                   lp_markers, iel_markers, skin_markers))

# the last row is called dummy and will be discarded  
markers_hit <- matrix(F, nrow = 1 + length(selected_features), ncol = ncol(all_markers), 
                      dimnames = list(c('dummy',selected_features), 
                                      colnames(all_markers)))

for(i in colnames(all_markers)){
  markers_hit[all_markers[,i],i] <- T
}
markers_hit <- as.data.frame(markers_hit)[-1,]

# top markers detected in tissues
head(rowSums(markers_hit[,1:10])[order(rowSums(markers_hit[,1:10]), 
                                       decreasing = T)])

# overlap_analysis --------------------------------------------------------
# compute intersection between marker sets and their significance (fisher test)
inters_mat = union_mat = fisher_mat = fisher_mat_bool = 
  matrix(nrow = ncol(all_markers), ncol = ncol(all_markers), 
         dimnames = list(colnames(all_markers), colnames(all_markers)))

for(row in rownames(inters_mat)){
  for(col in colnames(inters_mat)){
    inters_mat[row, col] <- sum(markers_hit[[row]] & markers_hit[[col]])
    union_mat[row, col] <- sum(markers_hit[[row]] | markers_hit[[col]])
    f_test <- fisher.test(markers_hit[,row], markers_hit[,col],
                          alternative = 'greater')
    p_value <- f_test$p.value
    fisher_mat[row,col] <- p_value 
    fisher_mat_bool[row,col] <- as.integer(p_value <0.01)
  }
}

# compute overlap between marker sets
overlap <- (inters_mat/union_mat)

# heatmap plots ------------------------------------------------------------
library(RColorBrewer)
cols = colorRampPalette(colors = c(brewer.pal(11,'Spectral'),'black'))(101)

overlap_hmap = pheatmap::pheatmap(overlap, 
                           color = rev(cols), 
                           breaks = (seq(0,1, length.out = 100)),
                           border_color = NA,
                           filename = './code/figures/6e_markers_overlap.pdf', 
                           fontsize = 5,annotation_legend = F,legend = F,
                           treeheight_row = 0, treeheight_col = 20,
                           width = 6.5/2.54, height = 7/2.54)

pheatmap::pheatmap(-log10(fisher_mat), 
                   color = c('black','khaki1', 'darkolivegreen1'),
                   number_color = 'black',
                   breaks = c(0, 2,5,20), border_color = NA,
                   cluster_rows = overlap_hmap[[1]], cluster_cols = overlap_hmap[[2]],
                   filename = './code/figures/6f_markers_overlap_pval.pdf', 
                   fontsize = 5,annotation_legend = F,legend = F,
                   treeheight_row = 0, treeheight_col = 20,
                   width = 6.5/2.54, height = 7/2.54)

