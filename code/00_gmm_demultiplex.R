library(ggplot2)
library(cowplot)
library(ggrastr)
library(mixtools)
library(mclust)
library(plyr)
library(EnvStats)
library(tidyverse)
set.seed(210201)

sample_files = c('./data/0_demultiplexing/00_s1_HTO_hashtag.count.tsv',
            './data/0_demultiplexing/00_s2_HTO_hashtag.count.tsv')

hto_demultiplexed = lapply(
  sample_files, 
  function(sample){
  hto <- read.table(sample, sep = "\t", row.names=1, header=T)
  
  # sum all hto counts for each barcode
  hto_sum <- data.frame(row.names=rownames(hto), total=rowSums(hto))
  
  # remove undetected hashtags
  hto <- hto[hto_sum>0, ]
  hto_sum <- hto_sum[hto_sum > 0,,drop=F]
  hto_freq <- hto / hto_sum$total
  
  # GMM model ---------------------------------------------------------------
  set.seed(1111)
  mcl.model <- Mclust(log10(hto_freq+1), G =  6) # model using six mixtures (# tissues +1)
  class.rename <- apply(mcl.model$parameters$mean, 1, which.max) # assign classes to indexes
  class = mapvalues(mcl.model$classification, from = class.rename, to = names(class.rename))
  unassigned_index <- intersect(unique(class), 1:6)
  
  # set stricter threshold for classification --------------------------------------------------
  max_prob <- apply(mcl.model$z[,-unassigned_index], MARGIN = 1, max)
  
  # set a threshold for assignment
  certain <- max_prob > 1-1e-4 
  # copy original classification
  classification2 <- mcl.model$classification
  # add z-score column
  z_score <- unlist(lapply(seq_along(mcl.model$classification),
                           FUN =  function(i) {mcl.model$z[i, mcl.model$classification[i]]}))
  # add uncertainity column
  uncertainty <- mcl.model$uncertainty
  # assign uncertain barcodes
  classification2[!certain] <- unassigned_index
  classification2[hto_sum$total < 200] <- unassigned_index
  classification2 <- mapvalues(classification2, from = class.rename, to = names(class.rename))
  
  class = mapvalues(classification2, from = class.rename, to = names(class.rename))
  # name non-named class as Unassigned
  class <- mapvalues(class, from = intersect(unique(class), 1:6), to = "Unassigned") 
  # extract sample number to add to rownames
  hto_stats <- data.frame(hto_freq, class=class, z_score = z_score, 
                                  uncertainty = uncertainty)
  sample_n = gsub('^.*0_demultiplexing/00_|_HTO_hashtag.*$','',sample)
  rownames(hto_stats) <- paste0(sample_n,'_',rownames(hto_stats),'-1')
  return(hto_stats)
})


write.csv(rbind(hto_demultiplexed[[1]],hto_demultiplexed[[2]]),
          file = './data/0_demultiplexing/01_HTO_gmm_demultiplexed.csv', quote = F)
