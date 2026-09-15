library(tidyverse)
library(readxl)
library(magrittr)

# read in subgraph membership
subgraph_membership <- read_xlsx("../../../5_phylogeny/3_generate_graph/subgraph_membership_final.xlsx")

# read in list of clusters
cluster_membership <- read_tsv("../../../2_protein_orf_processing/0.9_0.4_clusters.tsv",col_names = F)
colnames(cluster_membership) <- c("centroid","members")

# read in sequences with no host assignment to drop
no_host <- read_tsv("seqs_to_drop.txt",col_names=F)

write_drop_seqs <- function(target_subgraph){

  target_members <- subgraph_membership %>% 
    filter(Subgraph == target_subgraph)

  # drop any sequences belonging to clusters of any member in the blocked subgraph
  to_drop <- cluster_membership %>% 
    filter(centroid %in% target_members$Accession)
  
  write_lines(c(no_host$X1,to_drop$members),paste0("subgraph_",target_subgraph,"/all_drop_seqs.txt"))
  # query centroids only
  write_lines(paste0(target_members$Accession,"|Nucleocapsid"),paste0("subgraph_",target_subgraph,"/N_query_seqs.txt"))
  write_lines(paste0(target_members$Accession,"|Matrix"),paste0("subgraph_",target_subgraph,"/M_query_seqs.txt"))
  write_lines(paste0(target_members$Accession,"|Fusion"),paste0("subgraph_",target_subgraph,"/F_query_seqs.txt"))
  write_lines(paste0(target_members$Accession,"|Attachment"),paste0("subgraph_",target_subgraph,"/HN_query_seqs.txt"))
  write_lines(paste0(target_members$Accession,"|Polymerase"),paste0("subgraph_",target_subgraph,"/L_query_seqs.txt"))

}

lapply(unique(subgraph_membership$Subgraph),write_drop_seqs)
