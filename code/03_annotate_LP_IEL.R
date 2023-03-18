# includes ----------------------------------------------------------------
set.seed(210201)

source('./bin/useful_functions.R')
library(scater)
library(scran)
library(Seurat)

library(tidyverse)
library(ggthemes)
library(ggpubr)
library(patchwork)

# import ------------------------------------------------------------------

sce = readRDS('./data/1_preprocessing/01_all_tissues.rds')
rownames(sce) = uniquifyFeatureNames(ID = rowData(sce)[,'ID'], names = rowData(sce)[,'Symbol'])

# explorative plots -------------------------------------------------------
coldata_dr(sce, dimred = 'UMAP', color.column =  'tissue', size =1.5) + 
  scale_color_tableau('Tableau 20')

cols = c(tableau_color_pal('Tableau 20')(20),tableau_color_pal('Classic 10')(10) )
coldata_dr(sce[,sce$tissue %in% c( 'Colon','LP','IEL')], dimred = 'UMAP', color.column =  'global_louvain1.5', size =1.5) + 
  scale_color_manual(values = cols)

gene_dr(sce, dimred = 'UMAP', gene =  'Itgae', size =1.5) + scale_color_viridis_c()

# compute IEL, LP, tissue markers -----------------------------------------------------------------
tissue_markers_up = findMarkers(sce, groups = sce$tissue, 
                                direction = 'up', test.type = 'binom')
tissue_markers_dn = findMarkers(sce, groups = sce$tissue, 
                                direction = 'down', test.type = 'binom')

seu = as.Seurat(sce[,sce$tissue != 'IEL'])

int_markers = c('Cd244','Ccr9','Itgae', 'Asb2', # up in IEL
                 'Cd28' ,'Cd2','Ccr2')  # down in IEL

# IEL marker expressions in global dataset
Idents(seu) = factor(as.character(seu@meta.data[,'global_louvain1.5']))
DotPlot(seu, features = rev(int_markers), dot.scale = 7.5, 
        scale.by = 'size',cluster.idents = T) +
  RotatedAxis() + coord_flip()

# initial assignment --------------------------------------------------------------
# initially annotate as IEL that are not lp, in the refinement section this
# annotation will be used to select features to discern IEL from LP in SI and
# Colon
# 
iel_clusters = c(5,10,0,6,20,12,13) 

sce$tissue_hires = sce$tissue
sce$tissue_hires[sce$tissue == 'LP' & sce$global_louvain1.5 %in% iel_clusters] = 'LP_IEL'
sce$tissue_hires[sce$tissue == 'Colon' & sce$global_louvain1.5 %in% iel_clusters] = 'Colon_IEL'
sce$tissue_hires[sce$tissue == 'Colon' & (!sce$global_louvain1.5 %in% iel_clusters)] = 'Colon_LP'


# refinement --------------------------------------------------------------
## feature selection -----------------------------------------------------------------
# use initial assignment to obtain markers to discern IEL from LP

lp_iel = sce[,sce$tissue_hires %in% c('LP','IEL')]
tissue_markers = findMarkers(lp_iel, groups = lp_iel$tissue, 
                             test.type = 'binom')
genesel = rownames(tissue_markers$IEL[tissue_markers$IEL$FDR < 1e-10,])

## remove SI-IEL from SI-LP  dataset----------------------------------------------------------------------
lp = sce[,sce$tissue  == 'LP']
scran_clusters = quickCluster(lp, min.size = 100)
lp = computeSumFactors(lp, cluster = scran_clusters)
lp = logNormCounts(lp)

lp = runPCA(lp, subset_row = genesel, n_dimred = 20, name = 'pca_lp_iel')
lp = runUMAP(lp, subset_row = genesel,dimred = 'pca_lp_iel',
              name = 'umap_lp_iel')
lp = seurat_clusters(lp, dim.red = 'pca_lp_iel', 
                     resolution_param = 1.5, ndimred = 20, colname = 'lp_iel_clus')

# visualization of lp vs iel clusters
coldata_dr(lp, color.column = 'lp_iel_clus') + scale_color_tableau('Tableau 20')
coldata_dr(lp, color.column = 'lp_iel_clus', dimred = 'umap_lp_iel') + scale_color_tableau('Tableau 20')

seu = as.Seurat(lp)
Idents(seu) = factor(as.character(seu@meta.data[,'lp_iel_clus']))
DotPlot(seu, features = genesel[1:20],  dot.scale = 5,
        scale.by = 'size',cluster.idents = T) +
  RotatedAxis() + coord_flip()

# from marker expression, clusters 2,9,8 have an IEL-like profile
iel_clusters = c(2,9,8)
# 
lp$tissue_hires = lp$tissue
lp$tissue_hires[lp$tissue == 'LP' & lp$lp_iel_clus %in% iel_clusters] = 'LP_IEL'
lp$tissue_hires[lp$tissue == 'LP' & (!lp$lp_iel_clus %in% iel_clusters)] = 'LP'
# 
## split colon in LP and IEL --------------------------------------------------------------
colon = sce[,sce$tissue  == 'Colon']

scran_clusters = quickCluster(colon, min.size = 100)
colon = computeSumFactors(colon, cluster = scran_clusters)
colon = logNormCounts(colon)

colon = runPCA(colon, subset_row = genesel, n_dimred = 20, name = 'pca_lp_iel')
colon = runUMAP(colon, subset_row = genesel, dimred = 'pca_lp_iel', name = 'umap_lp_iel')
colon = seurat_clusters(colon, dim.red = 'pca_lp_iel', resolution_param = 1.5, ndimred = 20, colname = 'lp_iel_clus')

coldata_dr(colon, color.column = 'lp_iel_clus', size = 1.25, alpha =0.75) + 
  scale_color_tableau('Tableau 20') + 
  guides(colour = guide_legend(override.aes = list(size=5)))
coldata_dr(colon, color.column = 'lp_iel_clus',dimred = 'umap_lp_iel', size = 1.25) + 
  scale_color_tableau('Tableau 20')

seu = as.Seurat(colon)
Idents(seu) = factor(as.character(seu@meta.data[,'lp_iel_clus']))
DotPlot(seu, features = genesel[1:20],  dot.scale = 5,scale.by = 'size', cluster.idents = T,) +
  RotatedAxis() + coord_flip()

iel_clusters = c(3,0,7,11)

colon$tissue_hires = colon$tissue
colon$tissue_hires[colon$tissue == 'Colon' & colon$lp_iel_clus %in% iel_clusters] = 'Colon_IEL'
colon$tissue_hires[colon$tissue == 'Colon' & (!colon$lp_iel_clus %in% iel_clusters)] = 'Colon_LP'

# assign tissue --------------------------------------------------------------
colData(sce)[colnames(lp)[lp$tissue_hires == 'LP_IEL'], 'tissue_hires'] = 'LP_IEL'
colData(sce)[colnames(lp)[lp$tissue_hires == 'LP'], 'tissue_hires'] = 'LP'
colData(sce)[, 'lp_iel_clus'] = NA
colData(sce)[colnames(lp), 'lp_iel_clus'] = as.character(colData(lp)[,'lp_iel_clus'])

colData(sce)[colnames(colon)[colon$tissue_hires == 'Colon_IEL'], 'tissue_hires'] = 'Colon_IEL'
colData(sce)[colnames(colon)[colon$tissue_hires == 'Colon_LP'], 'tissue_hires'] = 'Colon_LP'

colData(sce)[colnames(colon), 'lp_iel_clus'] = as.character(colData(colon)[,'lp_iel_clus'])

coldata_dr(sce, color.column = 'tissue_hires',dimred = 'UMAP', size = 1.25) + 
  scale_color_tableau('Tableau 20')

# save results --------------------------------------------------------------------------------

saveRDS(sce, './data/1_preprocessing/02_all_tissues_hires.rds')


# {figure S6F} global UMAP by tissue (colon split) --------------------------------------------

#sce = readRDS('./data/1_preprocessing/02_all_tissues_hires.rds')

# plot all cells together and extract frames to use
# for each tissue plot
fig_S6F <- coldata_dr(sce, color.column =  'tissue_hires',size = 1) + 
  scale_color_tableau() + 
  guides(colour = guide_legend(override.aes = list(size=2)))

fig_S6F
ggsave('./code/figures/s6f_SI_lp_iel_annotation.pdf',
       plot = fig_S6F,
       width = 8, height = 4)

# {figure S6D} Small-intestine LP-IEL annotation -------------------------------------------------------------
# lp
seu = as.Seurat(lp)
Idents(seu) = factor(as.character(seu@meta.data[,'lp_iel_clus']))
lp_left = DimPlot(seu,reduction = 'umap_lp_iel',pt.size = 0.25) +
  scale_color_tableau('Tableau 20', direction =-1) + 
  scale_fill_tableau('Tableau 20', direction =-1) 

lp_right = DotPlot(seu, features = genesel[1:15], dot.scale = 2,scale.by = 'size',cluster.idents = T,) +
  RotatedAxis() + coord_flip() + guides(size =F, fill =F, col = F)

fig_S6D = lp_left + lp_right
ggsave('./code/figures/s6d_SI_lp_iel_annotation.pdf',
       plot = fig_S6D,
       width = 8, height = 4)

# {figure S6E} Small-intestine LP-IEL annotation -------------------------------------------------------------
# colon
seu = as.Seurat(colon)
Idents(seu) = factor(seu@meta.data[,'lp_iel_clus'])
colon_left = DimPlot(seu,reduction = 'umap_lp_iel',pt.size = 0.25) +
  scale_color_tableau('Tableau 20', direction =-1) + 
  scale_fill_tableau('Tableau 20', direction =-1)

colon_right = DotPlot(seu, features = genesel[1:15],  dot.scale = 2, scale.by = 'size',cluster.idents = T,) +
  RotatedAxis() + coord_flip() + guides(size =F, fill =F, col = F)

fig_S6E = colon_left + colon_right

ggsave('./code/figures/s6e_Colon_lp_iel_annotation.pdf',
       plot = fig_S6E,
       width = 8, height = 4)
