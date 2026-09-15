library(igraph)
library(readxl)
library(tidyverse)
library(writexl)
library(magrittr)
library(RColorBrewer)
library(rcartocolor)

### 
# this script calculates the mean vertex strength of blocked vertices
# this allows for plotting of relative increase in connectivity within the blocked subgraph as nodes are introduced
###

# read in graph structure
### read in phylogeny
p_phylo <- ape::read.nexus("../../5_phylogeny/2_beast_phylo/GTR_HMC_HIPSTR_tree.nxs")

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

# relabel tips with accessions instead of hosts
p_phylo$tip.label <- tip_hosts$ref_accessions

pat_distances <- ape::cophenetic.phylo(p_phylo)

### now generate a graph from this data
# it already functions as an adjacency matrix of a fully connected graph
g <- graph_from_adjacency_matrix(pat_distances, mode = "max", weighted = T)

# read in subgraph membership
subgraphs <- read_xlsx("../../5_phylogeny/deprecated/3_generate_graph/subgraph_membership_quant_15.xlsx")

# define a function to calculate a mean measure of 'closeness' across all blocked nodes
weighted_degree <- function(blocks){
  
  blocks <- blocks[colnames(blocks) != "subgraph7" & colnames(blocks) != "subgraph5"] # filter out minor subgraphs with membership n <= 5
  
  subgraph_degree <- vector(length = ncol(blocks))
  
  for(i in 1:ncol(blocks)){
    blocked_index <- blocks[,i][is.na(blocks[,i]) == F] # make a vector of the blocked indices
    
    blocked_degree <- vector(length = length(blocked_index))
    
    for(l in 1:length(blocked_index)){
      # loop over each blocked node, generate an minimum spanning tree using that node + all training nodes and calculate strength
      drop_nodes <- V(g)$name[blocked_index]
      
      g_out <- delete_vertices(g,drop_nodes[-l]) # drop all blocked nodes except for current target

      mg_out <- mst(g_out, weights = 1/E(g_out)$weight) # generate an mst from these nodes with reciprocal edge weights (> weight = weaker relation)
      # mg_out <- mst(g_out, weights = E(g_out)$weight) # generate an mst from these nodes
      
      # blocked_degree[l] <- strength(mg_out,vids = V(mg_out)$name[V(mg_out)$name == drop_nodes[l]]) # strength, or weighted vertex degree, measures the closeness to neighbouring vertices
      blocked_degree[l] <- harmonic_centrality(mg_out,vids = V(mg_out)$name[V(mg_out)$name == drop_nodes[l]]) # harmonic centrality measures closeness to all other vertices
    }
    
    subgraph_degree[i] <- mean(blocked_degree)
  }
  
  out_tbl <- tibble("mean_blocked_degree" = subgraph_degree,"subgraph" = colnames(blocks))
  
  return(out_tbl)
}

all_blocks <- list()

# read in blocked subgraph folds
blocked_folds <- read_csv("../../6_GNN/1_generate_validation_sets/blocked_cv_folds.csv")

all_blocks[1] <- list(blocked_folds)

setwd("../../6_GNN/1_generate_validation_sets")
filelist <- list.files(path = "../../6_GNN/1_generate_validation_sets", pattern = paste0("reintro_blocked_*"))

reintro_blocks <- lapply(filelist,read_csv)

all_blocks <- append(all_blocks,reintro_blocks) # final list including base blocked cv folds (no reintroductions) and reintro folds

all_out <- do.call(rbind,lapply(all_blocks,weighted_degree))

# now add identifiers of what blocking each result comes from
reintro_tbl <- do.call(rbind,str_extract_all(filelist,pattern = "\\d"))

reintro_values <- c(0,reintro_tbl[,1])
reintro_iteration <- c(0,reintro_tbl[,2])

all_out$num_introduced <- rep(reintro_values,each = nrow(all_out)/length(all_blocks))
all_out$iteration <- rep(reintro_iteration,each = nrow(all_out)/length(all_blocks))

all_out$subgraph <- factor(all_out$subgraph,levels = c("subgraph1","subgraph2","subgraph3","subgraph4",
                                                       "subgraph6","subgraph8","subgraph9","subgraph10"))

# basic plot to see how strength changes as we reintroduce nodes
all_out_plot <- all_out %>% 
  group_by(subgraph,num_introduced) %>% 
  summarise(mean_strength = mean(mean_blocked_degree))

palette <- carto_pal(10, "Safe")

ggplot(all_out_plot,aes(x = num_introduced,y = mean_strength,fill=subgraph))+
  geom_col() +
  theme_bw()+
  ylab("Mean centrality of blocked nodes") +
  xlab("# nodes reintroduced into blocked subgraph") +
  scale_fill_manual(values = palette) +
  facet_wrap(~subgraph,nrow = 4)

# ruh roh, strength basically doesn't correlate with reintroductions at the moment
# this is probably because in the fully connected graph, most nodes have many connections at short range
# when we drop nodes and draw a new minimum spanning tree, the distance to the next closest neighbour doesn't change much for most nodes
# mentally visualising the fully connected graph in 3d instead of the 2d minimum spanning tree makes this obvious
# it might be nice to actually have a 3d graph plot just to illustrate this point

# two options:
# think about other ways of defining subgraphs - the way we do it now does seem to correlate well with biological reality though (most subgraphs clearly comprise a particular viral clade)
# clustering by edge betweenness possibly sensible?

# think about other ways of defining 'connectedness' - difficult with a fully connected graph though (concepts like betweenness not really relevant)

# or just accept that this is true. the result makes sense thinking about the graph structure. At least you learned something

### SO 
# we can show this
# we can also plot the CATEGORICAL gain rather than quantitative gain in 'closeness' as we reintroduce
# in other words, plot the number of samples with the same host label in the training set and test sets
# this is probably where we are getting most performance gain from
# we could also plot the number in the same viral clades (e.g. how many jeilongviruses in the train set) and see if that correlates with performance


