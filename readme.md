Code to reproduce scRNA-seq results contained in 

> Sacirbegovic F, Günther M, Greco A, Zhao D, Wang X, Zhou M, Rosenberger S, Oberbarnscheidt MH, Held W, McNiff J, Jain D, Höfer T, Shlomchik WD. Graft-versus-host disease is locally maintained in target tissues by resident progenitor-like T cells. Immunity. 2023 Feb 14;56(2):369-385.e6. doi: 10.1016/j.immuni.2023.01.003. Epub 2023 Jan 30. PMID: 36720219.
> 

The repository contains two folders: `code` and `data`. 

- Scripts in `code/`  are numbered based on the order of execution. The working directory for running the scripts is the one containing both `code/` and `data/` . In order to run the analysis, please generate count matrices for GEX and HTO matrices using standard `Seurat::Read10X` or Bioconductor equivalent `DropletUtils::read10xCounts` pipeline.
- Data needs to be downloaded from GEO under the accession number **GSE207485**. In the `data/other/` folder, you’ll find a list of useful files for the analysis:
    - Information about the R version and relative libraries, in `data/other/renv.lock` . Instructions on how to use this file and recreate the R environment are [here](https://rstudio.github.io/renv/articles/renv.html)
    - `scGSEA_cycle` file contain cell cycle estimated as explained in [this paper](https://pubmed.ncbi.nlm.nih.gov/32783885/)
    - `dissociation_genes.txt` contains a list of genes associated with dissociation, they are used in the analysis to remove poor quality cells from the skin

### Scripts used to generate figures

| Figure | Content | script |
| --- | --- | --- |
| 6A | global UMAP by tissue (colon not split) | 02_preprocess_all_tissues |
| S6A | global UMAP by proliferation | 02_preprocess_all_tissues |
| S6B | global UMAP split by tissue (colon not split) | 02_preprocess_all_tissues |
| S6C | Heatmap: cluster vs tissue (colon not split) | 02_preprocess_all_tissues |
| 6B | Heatmap: tissue markers | 03.figure_tissue_markers |
| 6C | UMAP of tissues clustered | 05.figure_6c_tissues |
| 6D | Proliferation vs Tcf7 | 03.figure_tcf_proliferation |
| 6E | Markers overlap  | 06_markers_overlap |
| 6F | Markers overlap (p-value) | 06_markers_overlap |
| S6D | LPL-IEL annotation (Small intestine) - UMAP + clusters | 03_annotate_LP_IEL |
| S6E | LPL-IEL annotation (Colon) - UMAP + clusters | 03_annotate_LP_IEL |
| S6F | global UMAP by tissue (colon split) | 03_annotate_LP_IEL |
| S6G | Heatmap: cluster vs tissue (colon split) | 03.figure_tissue_markers |
