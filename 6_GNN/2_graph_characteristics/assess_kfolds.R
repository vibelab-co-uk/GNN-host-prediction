library(tidyverse)
library(igraph)
library(magrittr)

# read in kfolds
cv_folds_5 <- read_csv("../1_generate_validation_sets/cv_folds_5.csv")

# generate MST graph
# read in maximum likelihood phylogeny
p_phylo <- ape::read.nexus("../../5_phylogeny/2_beast_phylo/GTR_HMC_HPSTR_tree_50.nxs")

hosts <- read_csv("../../3_host_assignment/all_cluster_hosts.csv")

# create tibble of tip labels and bind in hosts
phylo_tips <- tibble(ref_accessions = p_phylo$tip.label)

tip_hosts <- left_join(phylo_tips, hosts)

tip_hosts %<>%
  mutate(Host_rank = case_when(is.na(Host_rank) == T ~ "Undetermined",
                               .default = Host_rank))

# label the phylogeny with hosts and drop tips with undetermined hosts
p_phylo$tip.label <- factor(tip_hosts$Host_rank,levels=c("Primates","Rodentia","Carnivora","Chiroptera","Artiodactyla","Perissodactyla","Scandentia","Pholidota","Eulipotyphla",
                                                         "Aves","Reptilia","Fish","Undetermined"))

p_phylo <- ape::drop.tip(p_phylo,"Undetermined")

tip_hosts %<>%
  filter(Host_rank != "Undetermined")


p_phylo$tip.label <- tip_hosts$ref_accessions
  
pat_distances <- ape::cophenetic.phylo(p_phylo)

g <- graph_from_adjacency_matrix(pat_distances, mode = "max", weighted = T)
mg <- igraph::mst(g, weights = igraph::E(g)$weight) # generate minimum spanning tree

# for each fold, calculate max and mean distance between blocked (test) nodes and train nodes

train_test_dist <- function(x){
  
  test_nodes <- cv_folds_5 %>% 
    select(colnames(cv_folds_5)[x]) %>% 
    drop_na()
  
  test_nodes <- data.frame(test_nodes)
  
  test_mask = vector(length = 212)
  
  test_mask[test_nodes[,1]] <- T
  
  train_mask <- !test_mask
  
  all_dist <- igraph::distances(graph = mg, v = V(mg)[test_mask], to = V(mg)[train_mask],weights = NA)
  
  query_nodes <- row.names(all_dist)
  
  all_dist <- as_tibble(all_dist)
  all_dist$query <- query_nodes  
  
  all_dist %<>%
    pivot_longer(cols = -query,names_to = "target",values_to = "path_length")
  
  all_dist %<>% 
    group_by(query) %>% 
    summarise(min = min(path_length))
  
  all_dist$fold <- x
  return(all_dist)
}

all_out <- do.call(rbind,lapply(1:5,train_test_dist))

table(all_out$min)
