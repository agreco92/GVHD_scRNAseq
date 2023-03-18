gene_dr <- function(sce, gene, dimred = 'UMAP', components = c(1,2), assay = 'logcounts', alpha=.8, size = 1){
  # ensembl <- rownames(sce)[rowData(sce)[,'Symbol'] == gene]
  # return(ensembl)
  df <- as.data.frame(reducedDim(sce, dimred)[,components])
  colnames(df) <- c('dr1','dr2')
  df[['gene']] <- assay(sce, assay)[gene,]
  # saturate at .99 quantile 
  df$gene[df$gene >= quantile(df$gene, .999)] = quantile(df$gene, .999)
  if(max(df$gene) == 0){return(ggplot(df) + geom_blank())}
  else{
    return(ggplot(df, aes(x = dr1, y = dr2)) + geom_point(aes(col = gene), alpha =alpha, size =size,stroke =0) + scale_color_viridis_c() +
               labs(x = paste0(dimred,'_',components[1]), y = paste0(dimred,'_',components[2]),colour = gene))
  }
}

coldata_dr <- function(sce, color.column, dimred = 'UMAP', components = c(1,2), alpha=.8, size = 1){
  #return(ensembl)
  df <- as.data.frame(reducedDim(sce, dimred)[,components])
  colnames(df) <- c('dr1','dr2')
  df[['color.column']] <- colData(sce)[,color.column]
  # saturate at .999 quantile 
  #df$color.column[df$color.column >= quantile(df$color.column, .999)] = quantile(df$color.column, .999)
  ggplot(df, aes(x = dr1, y = dr2)) + geom_point(aes(col = color.column), alpha =alpha, size =size,stroke =0) +
    labs(x = paste0(dimred,'_',components[1]), y = paste0(dimred,'_',components[2]), colour = color.column)
  }


seurat_clusters = function(sce, dim.red = 'PCA', resolution_param = 1.0, colname = 'cluster', ndimred = 30){
  sce.tmp <- sce
  rownames(sce.tmp) <- make.unique(rownames(sce.tmp))
  tmp <- as.Seurat(sce.tmp)
  tmp = FindNeighbors(tmp,reduction = dim.red, dims = 1:ndimred)
  tmp = FindClusters(tmp,resolution = resolution_param)
  colData(sce)[,colname] = tmp@active.ident
  return(sce)
}

ensembl_to_symbol <- function(sce, genes){
  return(rowData(sce)[genes,'Symbol'])}


convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  return(humanx)
}