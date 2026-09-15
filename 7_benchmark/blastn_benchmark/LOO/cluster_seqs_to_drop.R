library(tidyverse)
library(readxl)
library(magrittr)

# get a list of all directories (centroids)
directories <- list.files(path = ".",pattern = "\\.\\d$")

# read in list of clusters
cluster_membership <- read_tsv("../../../2_protein_orf_processing/0.9_0.4_clusters.tsv",col_names = F)
colnames(cluster_membership) <- c("centroid","members")

# read in sequences with no host assignment to drop
no_host <- read_tsv("seqs_to_drop.txt",col_names=F)

write_drop_seqs <- function(target_centroid){

  to_drop <- cluster_membership %>% 
    filter(centroid == target_centroid)

  write_lines(c(no_host$X1,to_drop$members),paste0(target_centroid,"/all_drop_seqs.txt"))
  write_lines(paste0(target_centroid,"|Nucleocapsid"),paste0(target_centroid,"/N_query_seqs.txt"))
  write_lines(paste0(target_centroid,"|Matrix"),paste0(target_centroid,"/M_query_seqs.txt"))
  write_lines(paste0(target_centroid,"|Fusion"),paste0(target_centroid,"/F_query_seqs.txt"))
  write_lines(paste0(target_centroid,"|Attachment"),paste0(target_centroid,"/HN_query_seqs.txt"))
  write_lines(paste0(target_centroid,"|Polymerase"),paste0(target_centroid,"/L_query_seqs.txt"))
  
}

lapply(directories,write_drop_seqs)
