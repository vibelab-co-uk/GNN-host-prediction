library(igraph)
library(ape)
library(tidyverse)
library(writexl)
library(magrittr)
library(RColorBrewer)
library(rcartocolor)

# read in maximum likelihood phylogeny
p_phylo <- ape::read.nexus("../2_beast_phylo/L_phylogeny/L_orfs_GTR_Strict_HIPSTR_summary.nxs")

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

# generate a graph from this data -----------------------------------------
#_______________________________________________________________________________
#_______________________________________________________________________________

# it already functions as the adjacency matrix of a fully connected graph
g <- graph_from_adjacency_matrix(pat_distances, mode = "max", weighted = T)

# set up a colour palette for plotting
# palette <- c(carto_pal(12, "Safe"),"grey1")
palette <- carto_pal(12, "Safe")

mycol = palette[p_phylo$tip.label]

# generate an mst graph ---------------------------------------------------
mg <- igraph::mst(g, weights = igraph::E(g)$weight) # generate minimum spanning tree

# take the reciprocal
rmg <- mg

E(rmg)$weight <- 1/E(rmg)$weight

write.csv(names(V(rmg)),"./graph_data/patristic_paramyxo_graph_node_names_L.csv")
write.csv(igraph::E(rmg)$weight,"./graph_data/patristic_paramyxo_graph_branch_weights_L.csv")
write_graph(rmg, "./graph_data/patristic_paraymxo_graph_L.edgelist", format = c("edgelist"))

# identify subgraphs --------------------
#_______________________________________________________________________________
#_______________________________________________________________________________

V(g)$name <- tip_hosts$ref_accessions

### manual subgraphs - find a threshold that preserves separate henipacluster with combined metaavulavirus cluster
quant <- quantile(igraph::E(g)$weight,probs = seq(0,1,0.01))

subgraph_edges <- c(igraph::E(g)[igraph::E(g)$weight<=2.14]) 
subgraph_edges <- subgraph_from_edges(g,subgraph_edges)

png(filename="./graph_phylo_plots/host_labels_subgraphs_final.png",width=1200,heigh=2400,units="px")
plot(subgraph_edges, vertex.size=3, vertex.label=NA,vertex.color=mycol)
dev.off()

# look at what sequences belong to which subclusters
subgraph_membership <- tibble("Accession" = names(components(subgraph_edges)$membership),"Subgraph" = components(subgraph_edges)$membership)

subgraph_membership %<>%
  arrange(Subgraph)

write_xlsx(subgraph_membership,"subgraph_membership_L.xlsx")
  
# subgraph 1 contains the clade at the bottom of the tree
# this is ortho- and para-rubulaviruses including mumps, infecting primarily bats and primates (humans)

# subgraph 2 contains the cluster above the birds
# this comprises reptilian and fish-infecting viruses and a number of mammalian respiroviruses (mostly human and livestock)

# subgraph 3 contains the primarily carnivora + artiodactyla cluster
# this comprises mammalian-infecting morbilliviruses

# subgraph 4 contains a cluster between respiroviruses and morbilliviruses
# this comprises mammalian-infecting narmoviruses (almost all rodent viruses + one bovine)

# subgraph 5 contains the cluster above the birds
# this comprises reptile-infecting viruses

# subgraph 6 contains the large top cluster
# this comprises mammalian-infecting jeilongviruses

# subgraph 7 is the orphan clade containing only 2 viruses, the Paraavulavirus wisconsinense
# these are avian infecting viruses on a very long branch from all the other avian viruses

# subgraph 8 is a small clade of pigeon-infecting metaavulaviruses

# subgraph 9 contains the upper part of the avian-infecting clade
# this comprises the avian orthoavulaviruses

# subgraph 11 contains some of the lower part of the avian-infecting clade
# this comprises avian metaavulaviruses

# subgraph 12 contains the henipaviruses

# subgraph 13 contains some of the lower part of the avian-infecting clade
# this comprises avian metaavulaviruses
##########
 