# choose resolution -------------------------------------------------------
set.seed(210201)

source('./bin/useful_functions.R')
library(scater)
library(scran)
library(Seurat)
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(patchwork)

tissues = c('Spleen','Skin_filtered', 'LP', 'IEL', 'Colon_LP','Colon_IEL')
clustering_res = c(0.5,0.3,0.4,0.4,0.4,0.75) %>% as.character()
names(clustering_res) = tissues

cluster_names = list(c('eff1','prog','cc1','cc2','eff2','ifn','dummy'), 
                     c('prog','eff2','eff1','cc2','dummy','dummy','dummy'),
                     c('prog','eff1','eff2','cc1','cc2','ifn','dummy'),
                     c('eff1','eff2','cc2','cc1','ifn','eff3','prog'),
                     c('eff1','eff2','eff3','prog','cc2','dummy','dummy'),
                     c('eff1','prog','eff2','cc2','dummy','dummy','dummy'))
names(cluster_names) = tissues

for(t in tissues){
  tissue_sce = readRDS(paste0('./data/2_tissue_preprocessing/',t,'.rds'))
  tissue_sce$cluster = 
    plyr::mapvalues(colData(tissue_sce)[,paste0('louvain',clustering_res[[t]])],
                    from =  c(0,1,2,3,4,5,6),
                    to = cluster_names[[t]]) %>% factor()
  
  colormap = c('#12a2a8', "#78a641", "#bcbd22", "#ffaa0e","#d63a3a","#ff7f0e","#c7519c")
  names(colormap) = c('prog', 'cc1','cc2','eff1','eff2','eff3','ifn')
  
  
  df <- data.frame(x1 = reducedDim(tissue_sce, 'UMAP')[,1], 
                   x2 = reducedDim(tissue_sce, 'UMAP')[,2], 
                   cluster = colData(tissue_sce)[,'cluster'])
  df_median <- df %>% group_by(cluster) %>% 
    summarise(x1 = median(x1, na.rm = T), x2 = median(x2, na.rm = T))
  
  dr_plot <- coldata_dr(tissue_sce, color.column = 'cluster', size = 1) + 
    geom_label(data = df_median, inherit.aes = F, 
               aes(x = x1, y = x2, label = cluster, col = cluster),alpha =.5)+
    scale_color_manual(values = colormap) +
    theme_few()
  
  # annotation heatmap ------------------------------------------------------
  seu = as.Seurat(tissue_sce)
  seu <- ScaleData(seu)
  Idents(seu) = factor(seu@meta.data[,'cluster'], 
                       levels = c('prog','eff1','eff2','eff3',
                                  'ifn','cc1','cc2'))
  
  int_markers <- c('Tcf7','Il7r','Tmsb10', # progenitor
                   'Gzma','Ccl5','Ccl4','Ccl3','Ifng', # effector
                   'Ifit1','Ifit3', #interferon
                   'Gapdh','Ppia','Ptma', # cell cycle 2
                   'Birc5','Ccna2') # cell cycle 1
  heatmap = DoHeatmap(seu, features = int_markers, group.colors = colormap) + 
    theme(axis.title = element_blank(),
          legend.position = 'none',
          legend.title = element_blank())
  
  heatmap
  fig = (dr_plot + heatmap) + plot_layout(widths = c(1.5,1))
  ggsave(paste0('./code/figures/6c/',t,'.pdf'), width = 8.75, height = 4.5)

}