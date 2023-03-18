set.seed(210201)
library(scater)
library(scran)
library(Seurat)
library(tidyverse)
library(ggthemes)
library(patchwork)
source('./code/useful_functions.R')
# includes ------------------------------------------------------------------------------------

library(DropletUtils)
library(AnnotationHub)

library(org.Mm.eg.db)
source('./code/useful_functions.R')

part1 <- read10xCounts('./data/1_preprocessing/00_count_matrix/part1/filtered_feature_bc_matrix/',version = '3')
part2 <- read10xCounts('./data/1_preprocessing/00_count_matrix/part2/filtered_feature_bc_matrix/',version = '3')
demulti_gmm <- read.csv('./data/0_demultiplexing/01_HTO_gmm_demultiplexed.csv', row.names = 1)
demulti_bioc = read.csv(
  './data/0_demultiplexing/02_HTO_dropletutils_demultiplexed.csv', row.names = 1)

# join SingleCellExperiments and add Sample information
part1$Sample = 'Sample_1'
part2$Sample = 'Sample_2'
sce <- cbind(part1, part2)
colnames(sce) <- paste0(sce$Sample, '_', sce$Barcode)

# add GMM demultiplexing info -----------------------------------------------------------------
sce$tissue = 'Unassigned'
demulti_gmm_columns <- c('tissue','demultiplex_zscore','demultiplex_uncertainty')

colData(sce)[
  gsub('^s','Sample_', rownames(demulti_gmm)),demulti_gmm_columns] = 
  demulti_gmm[,c('class','z_score','uncertainty')]

# character to numerical conversion for relevant metadata
sce$demultiplex_zscore = as.numeric(sce$demultiplex_zscore)
sce$demultiplex_uncertainty = as.numeric(sce$demultiplex_uncertainty)

# add HashedDrops demultiplexing ---------------------------------------------------------
rownames_demulti = gsub('^Sample_','s', colnames(sce)) # rename cells 

colData(sce)[,c('Second', 'Confident')] = 
  demulti_bioc[rownames_demulti, c('Second', 'Confident')]

# quality control ---------------------------------------------------------
# remove cells assigned with low confidence
sce <- sce[,sce$tissue != 'Unassigned' & sce$Confident == T]

is.mito <- grep('^mt.*', rowData(sce)[,'Symbol']) #  mitochondrial transcripts
stats <- perCellQCMetrics(sce, subsets=list(Mt=is.mito))

# detect outliers
# sample2 from skin and spleen have a large number of bad cells,
# thus are not used for outlier computation 

discard_libsize <- isOutlier(stats$sum, type = 'lower', log = T, batch = sce$tissue, 
                             subset = !paste0(sce$tissue,'_',sce$Sample) %in% c('Skin_Sample_2','Spleen_Sample_2'), 
                             share_medians = T,share_mads = T)
discard_features <- isOutlier(stats$detected, type = 'lower', log = T, batch = sce$tissue, 
                              subset = !paste0(sce$tissue,'_',sce$Sample) %in% c('Skin_Sample_2','Spleen_Sample_2'), 
                              share_medians = T, share_mads = T)
discard_mito <- isOutlier(stats$subsets_Mt_percent, type = 'higher', log = F, batch = sce$tissue,
                          subset = !paste0(sce$tissue,'_',sce$Sample) %in% c('Skin_Sample_2','Spleen_Sample_2'), 
                          share_medians = T, share_mads = T)

# discard outliers and cells of low quality
discard_libsize <- discard_libsize | stats$sum < 1000
discard_features <- discard_features | stats$detected < 200
discard_mito <- discard_mito | stats$subsets_Mt_percent > 10

discard_df <- data.frame(libsize = discard_libsize, features = discard_features, mitoc = discard_mito, tissue = sce$tissue)
discard_df$discard <- sign(rowSums(discard_df[,c(1,2,3)]))

# write number of discarded cells
# table(discard_df$discard, discard_df$tissue) %>%
#   as.data.frame.array() %>% write.csv('./data/1_preprocessing/poor_quality_cells.csv', quote = F)

# add QC metadata
colData(sce)[,colnames(stats)] <- stats
sce$discard <- as.logical(discard_df$discard)

qc_plot <- gridExtra::grid.arrange(
  plotColData(sce, x="Sample", y="sum", other_fields = 'tissue', colour_by="discard", point_size =.5) +   facet_wrap(~tissue, nrow =2) + scale_y_log10() + ggtitle("Total count"),
  plotColData(sce, x="Sample", y="detected", other_fields = 'tissue', colour_by="discard" ,point_size =.5)  +   facet_wrap(~tissue, nrow =2) + scale_y_log10() + ggtitle("Detected features"),
  plotColData(sce, x="Sample", y="subsets_Mt_percent", other_fields = 'tissue', colour_by="discard", point_size =.5)  +   facet_wrap(~tissue, nrow =2) + ggtitle("Mt percent"),
  ncol=1)

sce <- sce[,!discard_df$discard]


# pre-processig -------------------------------------------------------------------------------

# compute sct
seu <- Seurat::CreateSeuratObject(counts = counts(sce))
seu <- SCTransform(seu, return.only.var.genes = T, verbose = T)

# dimensionality reduction
## PCA
sce_sct <- SingleCellExperiment(
  assays = list(counts = seu@assays$SCT@counts[seu@assays$SCT@var.features,],
                logcounts = seu@assays$SCT@scale.data[seu@assays$SCT@var.features,]))

sce_sct <- runPCA(sce_sct, subset_row = seu@assays$SCT@var.features)
## UMAP
reducedDim(sce, 'PCA') <- reducedDim(sce_sct)
sce_sct <- runUMAP(sce_sct,dimred = 'PCA', ncomponents =2, n_dimred =20)

# clustering
for(res in c(0.4,0.5,1)){
  sce_sct <- seurat_clusters(sce_sct,dim.red = 'PCA', ndimred = 20, 
                             resolution_param = res, colname = paste0('global_louvain', res))
}

# remove outlier cluster -----------------------------------------------------------------------
coldata_dr(sce_sct, 'global_louvain0.5', size =1) + 
  scale_color_tableau('Tableau 20')

sce <- sce[,sce_sct$global_louvain0.5 != 12]

# repeat preproccessing without outliers -----------------------------------------------------------------------
# compute sct
seu <- Seurat::CreateSeuratObject(counts = counts(sce))
seu <- SCTransform(seu, return.only.var.genes = T, verbose = T)

## store sct results ---------------------------------------------------------------------------
saveRDS(seu@assays$SCT@scale.data, './data/1_preprocessing/01_sctransform_residuals.rds')
saveRDS(seu@assays$SCT@var.features, './data/1_preprocessing/01_sctransform_hvgs.rds')
# dimensionality reduction
## PCA
sce_sct <- SingleCellExperiment(
  assays = list(logcounts = seu@assays$SCT@scale.data))

sce_sct <- runPCA(sce_sct, subset_row = seu@assays$SCT@var.features)
reducedDim(sce, 'PCA') <- reducedDim(sce_sct) # import in main object
## UMAP
sce <- runUMAP(sce,dimred = 'PCA', ncomponents = 2, n_dimred =20)

# Seurat's clustering function needs logcounts (although clustering is performed
# using PCA), counts will be normalized using the pooling method later
# in the analysis
sce <- logNormCounts(sce) 

for(res in c(0.25,0.4,0.5,0.75,1,1.5,2)){
  sce <- seurat_clusters(sce, dim.red = 'PCA', ndimred = 20, resolution_param = res, colname = paste0('global_louvain', res))
}

# {figure 6A} - global umap ---------------------------------------------------------------------
# export figure
p <- coldata_dr(sce, color.column =  'tissue',size = 1) + theme_few() + 
  scale_color_tableau() + # color scale
  guides(colour = guide_legend(override.aes = list(size=2))) + # big legend
  theme(axis.text = element_blank())

ggsave('./code/figures/1a_global_umap.pdf', p, width = 6, height = 4.5)

# scran normalization -----------------------------------------------------
scran_clusters <- quickCluster(sce, min.size = 500)
sce <- computeSumFactors(sce, cluster = scran_clusters)
sce <- logNormCounts(sce)

# estimate cell cycle -----------------------------------------------------

library(stringr)
seu <- as.Seurat(sce)
s.genes <- str_to_title(cc.genes.updated.2019$s.genes)
s.genes <- rownames(rowData(sce)[rowData(sce)$Symbol %in% s.genes,])

g2m.genes <- str_to_title(cc.genes.updated.2019$g2m.genes)
g2m.genes <- rownames(rowData(sce)[rowData(sce)$Symbol %in% g2m.genes,])

seu <- CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
colData(sce)[,c('seurat_G2M','seurat_S', 'seurat_ccphase')] <- seu@meta.data[,c('G2M.Score','S.Score','Phase')]

s_g2m.genes <- c(s.genes, g2m.genes)
seu <- AddModuleScore(seu, features = list(s_g2m.genes),name = 'cycle.Score')
colData(sce)[,'seurat_nonG1'] <- seu@meta.data$cycle.Score1

# {figure S6A} - proliferation -----------------------------------------------------------------------
p <- coldata_dr(sce, color.column =  'seurat_nonG1',size = 1.25) + theme_few() + 
  #scale_color_tableau() + # color scale
  guides(colour = guide_legend(override.aes = list(size=2))) + # big legend
  theme(axis.text = element_blank(),
        axis.title = element_text(size = 5),
        legend.text = element_text(size =5),
        legend.title = element_text(size =7),
        legend.key.size = unit(0, "pt"),
        legend.margin=margin(0,0,0,0, unit='pt'))

ggsave('./code/figures/6a_global_umap_proliferation.pdf', 
       p, width = 6, height = 4.5)


# {figure S6B} - global umap split by tissue --------------------------------------------------
sce = readRDS('./data/1_preprocessing/01_all_tissues.rds')
colors = tableau_color_pal('Tableau 10')(5)
names(colors) <- names(table(sce$tissue))

# plot all cells together and extract frames to use
# for each tissue plot
p <- coldata_dr(sce, color.column =  'tissue',size = 0.5) + 
  scale_color_manual(values = colors) + # color scale
  guides(colour = guide_legend(override.aes = list(size=2)))

xlims <- ggplot_build(p)$layout$panel_scales_x[[1]]$range$range
ylims <- ggplot_build(p)$layout$panel_scales_y[[1]]$range$range

theme_update(axis.title = element_blank())

p1 <- coldata_dr(sce[,sce$tissue == 'Spleen'],size = .75, color.column = 'tissue')  + scale_color_manual(values = colors) + guides(col =F) + xlim(xlims) + ylim(ylims)
p2 <- coldata_dr(sce[,sce$tissue == 'Skin'], size = .75,  color.column = 'tissue')  + scale_color_manual(values = colors) + guides(col =F) + xlim(xlims) + ylim(ylims)
p3 <- coldata_dr(sce[,sce$tissue == 'LP'], size = .75, color.column = 'tissue')  + scale_color_manual(values = colors) + guides(col =F) + xlim(xlims) + ylim(ylims)
p4 <- coldata_dr(sce[,sce$tissue == 'IEL'], size = .75, color.column = 'tissue') + scale_color_manual(values = colors) + guides(col =F) + xlim(xlims) + ylim(ylims)
p5 <- coldata_dr(sce[,sce$tissue == 'Colon'], size = .75, color.column = 'tissue')  + scale_color_manual(values = colors)+ guides(col =F) + xlim(xlims) + ylim(ylims)

p <- p1 + p2 + p3+ p4+ p5 + 
  plot_layout(ncol=3, guides = 'collect')
ggsave('./code/figures/s1b_global_umap_split_by_tissue.pdf', 
       p, width = 6, height = 4.5)

# {figure S6C} - tissue vs cluster --------------------------------------------------

cluster_tissue_confusion <- table(sce$global_louvain1, sce$tissue)
cluster_tissue_confusion <- cluster_tissue_confusion/rowSums(cluster_tissue_confusion)

cols <- viridis::plasma(101)

pheatmap::pheatmap(cluster_tissue_confusion, display_numbers = T, border_color = NA, 
                   color = cols, legend = F,
                   filename = './code/figures/s6c_tissue_vs_cluster.pdf',
                   width = 5, height = 4)


# save SingleCellExperiment object ------------------------------------------------------------

saveRDS(sce, './data/1_preprocessing/01_all_tissues.rds')
