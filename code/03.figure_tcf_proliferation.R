set.seed(210201)

library(scater)
library(scran)
library(ggthemes)
library(ggrepel)
library(patchwork)
library(tidyverse)
source('./bin/useful_functions.R')

sce <- readRDS('./data/1_preprocessing/02_all_tissues_hires.rds')

tissue_files = paste0('./data/2_tissue_preprocessing/',
                      c('Colon_IEL','Colon_LP','IEL','LP','Skin_filtered','Spleen'),
                      '.rds')

all_cells = lapply(
  tissue_files, 
  FUN = function(x){
    return(colnames(readRDS(x)))
  }) %>% Reduce(f = union, x =.)

sce <- sce[,all_cells]


# add more precise cycling estimate computed using scGSEA ----------------------------------------------
# substitute with cell cycle estimate from Seurat to obtain an equivalent plot

scgsea1 <- readRDS('./data/other/scGSEA_cycle_sample1.rds')
scgsea2 <- readRDS('./data/other/scGSEA_cycle_sample2.rds')


new_cols <- c('scgsea_G2M', 'scgsea_S', 'scgsea_ccphase','scgsea_nonG1')
tmp <- left_join(as.data.frame(colData(sce)[sce$Sample == 'Sample_1','Barcode', drop =F]), 
                 scgsea1 %>% rownames_to_column('Barcode'), by = 'Barcode')
colData(sce)[sce$Sample == 'Sample_1', new_cols] = tmp[,colnames(scgsea1)]

tmp <- left_join(as.data.frame(colData(sce)[sce$Sample == 'Sample_2','Barcode',drop =F]), 
                 scgsea2 %>% rownames_to_column('Barcode'), by = 'Barcode')
colData(sce)[sce$Sample == 'Sample_2', new_cols] = tmp[,colnames(scgsea2)]

sce <- sce[,!is.na(sce$scgsea_ccphase)]

cycle_vs_tcf7 <- data.frame(tcf7 = logcounts(sce)['Tcf7', ], 
                  tissue = sce$tissue_hires, 
                  phase = sce$scgsea_ccphase) %>% 
  group_by(tissue) %>% 
  summarise(tcf7_pos = sum(tcf7>0)/n(), 
            cycle = 1 - sum(phase == 'G1')/n(),
            cycle_tcf7 = sum(tcf7 > 0 & phase!='G1')/sum(tcf7>0))

fig_6D <- ggplot(cycle_vs_tcf7, aes(x = 100* tcf7_pos, y = 100* cycle_tcf7,  label = tissue)) + 
  geom_point() + 
  geom_label_repel(size = 5*0.352777778) + # points to mm (geom_label uses those to indicate size))
  ylab('Cycling Tcf7+ (%)') + xlab('Tcf7+ in tissue (%)') + 
  theme_few() + 
  theme(legend.position = 'none', text = element_text(size = 5), axis.text = element_text(size = 7),
        plot.margin = unit(c(0,0,0,0),'mm'))

p
ggsave('./code/figures/6d_tcf7_vs_proliferation.pdf',fig_6D, 
       width = 3.5, height = 3.5, units = 'cm')
