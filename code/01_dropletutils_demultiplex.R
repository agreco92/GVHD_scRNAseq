library(DropletUtils)
library(tidyverse)
library(scater)
set.seed(210201)

sample_files = c('./data/0_demultiplexing/00_hto_unfiltered/00_s1_HTO_hashtag.count.tsv',
                 './data/0_demultiplexing/00_hto_unfiltered/00_s2_HTO_hashtag.count.tsv')

hto_demultiplexed = lapply(
  sample_files,
  function(sample){
    hto =  read.table(sample, sep = "\t", row.names=1, header=T)
    
    # remove NA col, transpose
    hto = t(hto[,-ncol(hto)])
    hto_calls =  emptyDrops(hto, lower=200)
    
    is.cell = which(hto_calls$FDR <= 0.001)
    hto_stats = hashedDrops(hto,ambient=metadata(hto_calls)$ambient)
    
    # extract sample number to add back to barcode later
    sample_n = gsub('^.*00_hto_unfiltered/00_|_HTO_hashtag.*$','',sample)
    rownames(hto_stats) <- paste0(sample_n,'_',rownames(hto_stats),'-1')
    
    hto_stats$Best <- rownames(hto)[hto_stats$Best]
    hto_stats$Second <- rownames(hto)[hto_stats$Second]
    return(hto_stats)
  }) 

write.csv(rbind(hto_demultiplexed[[1]],hto_demultiplexed[[2]]),
          file = './data/0_demultiplexing/02_HTO_dropletutils_demultiplexed.csv', quote = F)