library(tidyverse)
library(readxl)
library(magrittr)

`%nin%` <- Negate(`%in%`)

# read in subgraph membership
subgraph_membership <- read_xlsx("../../../5_phylogeny/3_generate_graph/subgraph_membership_final.xlsx")

# read in list of clusters
cluster_membership <- read_tsv("../../../2_protein_orf_processing/0.9_0.4_clusters.tsv",col_names = F)
colnames(cluster_membership) <- c("centroid","members")

# read in graph-ordered node names
node_names <- read_csv("../../../5_phylogeny/3_generate_graph/graph_data/patristic_paramyxo_graph_node_names.csv")
colnames(node_names) <- c("Index","Accession")

# read in sequences with no host assignment to drop
no_host <- read_tsv("seqs_to_drop.txt",col_names=F)

write_drop_seqs <- function(reintro){

  n <- str_extract_all(reintro,"\\d")[[1]][1]
  iter <- str_extract_all(reintro,"\\d")[[1]][2]
  
  # read in reintro seqs
  train_indices <- read_csv(paste0("../../../6_GNN/1_generate_validation_sets/reintro_blocked_cv_folds_n",n,"_iteration",iter,".csv"))

  for(i in 1:length(colnames(train_indices))){
    
    train_seqs <- node_names %>% 
      filter(Index %in% as.vector(unname(train_indices[,i]))[[1]])
    
    # drop any sequences belonging to clusters not in the train set
    to_drop <- cluster_membership %>% 
      filter(centroid %in% train_seqs$Accession)
    
    write_lines(c(no_host$X1,to_drop$members),paste0(reintro,"/subgraph_",str_extract(colnames(train_indices[,i]),"\\d+"),"/all_drop_seqs.txt"))
    # query centroids only
    write_lines(paste0(train_seqs$Accession,"|Nucleocapsid"),paste0(reintro,"/subgraph_",str_extract(colnames(train_indices[,i]),"\\d+"),"/N_query_seqs.txt"))
    write_lines(paste0(train_seqs$Accession,"|Matrix"),paste0(reintro,"/subgraph_",str_extract(colnames(train_indices[,i]),"\\d+"),"/M_query_seqs.txt"))
    write_lines(paste0(train_seqs$Accession,"|Fusion"),paste0(reintro,"/subgraph_",str_extract(colnames(train_indices[,i]),"\\d+"),"/F_query_seqs.txt"))
    write_lines(paste0(train_seqs$Accession,"|Attachment"),paste0(reintro,"/subgraph_",str_extract(colnames(train_indices[,i]),"\\d+"),"/HN_query_seqs.txt"))
    write_lines(paste0(train_seqs$Accession,"|Polymerase"),paste0(reintro,"/subgraph_",str_extract(colnames(train_indices[,i]),"\\d+"),"/L_query_seqs.txt"))
  }
}

dirlist <- list.files(path = ".",pattern = "reintro_blocked_cv_n*")

lapply(dirlist,write_drop_seqs)
