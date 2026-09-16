library(igraph)
library(ape)
library(tidyverse)
library(writexl)
library(magrittr)
library(RColorBrewer)
library(rcartocolor)

# read in maximum likelihood phylogeny
p_phylo <- ape::read.nexus("../2_beast_phylo/GTR_HMC_HPSTR_tree_50.nxs")

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

pat_distances <- ape::cophenetic.phylo(p_phylo)

# relabel tips with accessions instead of hosts
p_phylo$tip.label <- tip_hosts$ref_accessions

write.tree(p_phylo,"GTR_HMC_HPSTR_tree_50_pruned.tree")

# generate a graph from this data -----------------------------------------
#_______________________________________________________________________________
#_______________________________________________________________________________

# it already functions as the adjacency matrix of a fully connected graph
g <- graph_from_adjacency_matrix(pat_distances, mode = "max", weighted = T)

# set up a colour palette for plotting
# palette <- c(carto_pal(12, "Safe"),"grey1")
palette <- carto_pal(12, "Safe")

mycol = palette[factor(tip_hosts$Host_rank,levels=c("Primates","Rodentia","Carnivora","Chiroptera","Artiodactyla","Perissodactyla","Scandentia","Pholidota","Eulipotyphla",
                                                    "Aves","Reptilia","Fish"))]

# read in graph layout
manual_layout <- read_csv("./graph_phylo_plots/final_manual_layout.csv")
manual_layout %<>% select(x_coord,y_coord)
manual_layout <- as.matrix(manual_layout)

# write a fully connected graph with exponentially scaled edgeweights --------
fg <- g
E(fg)$weight <- 1/(exp(2*E(fg)$weight)) # note we take the reciprocal to weight message passing correctly

write.csv(names(V(fg)),"./graph_data/patristic_paramyxo_expfull_graph_node_names.csv")
write.csv(igraph::E(fg)$weight,"./graph_data/patristic_paramyxo_expfull_graph_branch_weights.csv")
write_graph(fg, "./graph_data/patristic_paraymxo_expfull_graph.edgelist", format = c("edgelist"))

# generate an mst graph ---------------------------------------------------
mg <- igraph::mst(g, weights = igraph::E(g)$weight) # generate minimum spanning tree

# take the reciprocal
rmg <- mg

E(rmg)$weight <- 1/E(rmg)$weight

write.csv(names(V(rmg)),"./graph_data/patristic_paramyxo_graph_node_names.csv")
write.csv(igraph::E(rmg)$weight,"./graph_data/patristic_paramyxo_graph_branch_weights.csv")
write_graph(rmg, "./graph_data/patristic_paraymxo_graph.edgelist", format = c("edgelist"))

plot(mg, vertex.size=3, vertex.label=NA,vertex.color=mycol,layout = manual_layout)

# generate an mst graph with uniform edge weights -------------------------
ug <- mg
E(ug)$weight <- 1

write.csv(names(V(ug)),"./graph_data/patristic_paramyxo_unweighted_graph_node_names.csv")
write.csv(igraph::E(ug)$weight,"./graph_data/patristic_paramyxo_unweighted_graph_branch_weights.csv")
write_graph(ug, "./graph_data/patristic_paraymxo_unweighted_graph.edgelist", format = c("edgelist"))

# generate an mst graph with exponentially scaled edge lengths ------------
# exp_plot <- pivot_longer(tibble("linear_weights"=E(mg)$weight,"exponential_weights"=exp(E(mg)$weight)),cols=everything(),names_to = "transformation",values_to = "weights")
# 
# ggplot(exp_plot,aes(y = weights,x = transformation))+
#   geom_point()

eg <- mg
E(eg)$weight <- 1/(exp(E(mg)$weight))

write.csv(names(V(eg)),"./graph_data/patristic_paramyxo_exp_graph_node_names.csv")
write.csv(igraph::E(eg)$weight,"./graph_data/patristic_paramyxo_exp_graph_branch_weights.csv")
write_graph(eg,"./graph_data/patristic_paraymxo_exp_graph.edgelist", format = c("edgelist"))

# generate alternative graphs with greater connectivity within clusters --------
# keep the edges needed for the mst and drop all other edges longer than an arbitrary value
# look at some summary stats to select sensible value
edge_quant <- quantile(E(g)$weight,probs = seq(0,1,0.05))

additive_edges_5 <- c(E(g)[E(mg)],E(g)[E(g)$weight<=edge_quant[2]]) # keeping edges within region of 5th or 10th percentile looks biologically sensible
additive_edges_10 <- c(E(g)[E(mg)],E(g)[E(g)$weight<=edge_quant[3]])

ag_5 <- subgraph_from_edges(g,additive_edges_5)
ag_10 <- subgraph_from_edges(g,additive_edges_10)

png(filename="./graph_phylo_plots/host_labels_additive_graph_5.png",width=1200,heigh=2400,units="px")
plot(ag_5, vertex.size=3, vertex.label=NA,vertex.color=mycol)
dev.off()

png(filename="./graph_phylo_plots/host_labels_additive_graph_5_named.png",width=1200,heigh=2400,units="px")
plot(ag_5, vertex.size=3, vertex.color=mycol)
dev.off()

png(filename="./graph_phylo_plots/host_labels_additive_graph_10.png",width=1200,heigh=2400,units="px")
plot(ag_10, vertex.size=3, vertex.label=NA,vertex.color=mycol)
dev.off()

png(filename="./graph_phylo_plots/host_labels_additive_graph_10_named.png",width=1200,heigh=2400,units="px")
plot(ag_10, vertex.size=3, vertex.color=mycol)
dev.off()

rag_5 <- ag_5

E(rag_5)$weight <- 1/E(ag_5)$weight

write.csv(names(V(rag_5)),"./graph_data/patristic_paramyxo_additive_graph_5_node_names.csv")
write.csv(igraph::E(rag_5)$weight,"./graph_data/patristic_paramyxo_additive_graph_5_branch_weights.csv")
write_graph(rag_5, "./graph_data/patristic_paraymxo_additive_graph_5.edgelist", format = c("edgelist"))

rag_10 <- ag_10

E(rag_10)$weight <- 1/E(ag_10)$weight

write.csv(names(V(rag_10)),"./graph_data/patristic_paramyxo_additive_graph_10_node_names.csv")
write.csv(igraph::E(rag_10)$weight,"./graph_data/patristic_paramyxo_additive_graph_10_branch_weights.csv")
write_graph(rag_10, "./graph_data/patristic_paraymxo_additive_graph_10.edgelist", format = c("edgelist"))

# generate a more connected graph with exponential edge weights -----------
reag_5 <- ag_5
E(reag_5)$weight <- 1/(exp(E(ag_5)$weight))

write.csv(names(V(reag_5)),"./graph_data/patristic_paramyxo_exp_additive_graph_5_node_names.csv")
write.csv(igraph::E(reag_5)$weight,"./graph_data/patristic_paramyxo_exp_additive_graph_5_branch_weights.csv")
write_graph(reag_5, "./graph_data/patristic_paraymxo_exp_additive_graph_5.edgelist", format = c("edgelist"))


# generate a set of random spanning trees with exponential edge weights --------
for(i in 1:5){
  rand_st <- sample_spanning_tree(g, vid = 0)
  sg <- subgraph_from_edges(g,rand_st)
  #plot(sg, vertex.size=3, vertex.color=mycol)
  E(sg)$weight <- 1/(exp(E(sg)$weight))
  
  write.csv(names(V(sg)),paste0("./graph_data/random_spanning_trees/patristic_paramyxo_spanning_",i,"_graph_node_names.csv"))
  write.csv(igraph::E(sg)$weight,paste0("./graph_data/random_spanning_trees/patristic_paramyxo_spanning_",i,"_graph_branch_weights.csv"))
  write_graph(sg, paste0("./graph_data/random_spanning_trees/patristic_paraymxo_spanning_",i,"_graph.edgelist"), format = c("edgelist"))
}

# identify subgraphs --------------------
#_______________________________________________________________________________
#_______________________________________________________________________________

quant <- quantile(igraph::E(g)$weight,probs = seq(0,1,0.05))

V(g)$name <- tip_hosts$ref_accessions

### 15% subgraphs
subgraph_edges <- c(igraph::E(g)[igraph::E(g)$weight<=quant[4]]) 
subgraph_edges <- subgraph_from_edges(g,subgraph_edges)

png(filename="./graph_phylo_plots/host_labels_subgraphs_quant_15.png",width=1200,heigh=2400,units="px")
plot(subgraph_edges, vertex.size=3, vertex.label=NA,vertex.color=mycol)
dev.off()

# look at what sequences belong to which subclusters
subgraph_membership <- tibble("Accession" = names(components(subgraph_edges)$membership),"Subgraph" = components(subgraph_edges)$membership)

subgraph_membership %<>%
  arrange(Subgraph)

write_xlsx(subgraph_membership,"subgraph_membership_quant_15.xlsx")

# subgraph 1 contains the clade at the bottom of the tree
# this is ortho- and para-rubulaviruses including mumps, infecting primarily bats and primates (humans)

# subgraph 2 contains the cluster above the birds
# this comprises reptilian and fish-infecting viruses and a number of mammalian respiroviruses (mostly human and livestock)

# subgraph 3 contains the primarily carnivora + artiodactyla cluster
# this comprises mammalian-infecting morbilliviruses

# subgraph 4 contains the large top cluster
# this comprises mammalian-infecting jeilongviruses

# subgraph 5 contains the cluster above the birds
# this comprises reptile-infecting viruses

# subgraph 6 contains a cluster between respiroviruses and morbilliviruses
# this comprises mammalian-infecting narmoviruses (almost all rodent viruses + one bovine)

# subgraph 7 is the orphan clade containing only 2 viruses, the Paraavulavirus wisconsinense
# these are avian infecting viruses on a very long branch from all the other avian viruses

# subgraph 8 is a small clade of divergent paraavulavirus taiwanense
# these are pigeon infecting viruses on a very long branch from all the other taiwanense

# subgraph 9 contains the upper part of the avian-infecting clade
# this comprises the avian orthoavulaviruses

# subgraph 10 contains the lower part of the avian-infecting clade
# this comprises the avian metaavulaviruses

# subgraph 11 contains the mammalian cluster in between the morbilliviruses and narmoviruses
# this comprises the henipaviruses

### 20% subgraphs
subgraph_edges <- c(igraph::E(g)[igraph::E(g)$weight<=quant[5]]) 
subgraph_edges <- subgraph_from_edges(g,subgraph_edges)

png(filename="./graph_phylo_plots/host_labels_subgraphs_quant_20.png",width=1200,heigh=2400,units="px")
plot(subgraph_edges, vertex.size=3, vertex.label=NA,vertex.color=mycol)
dev.off()

# look at what sequences belong to which subclusters
subgraph_membership <- tibble("Accession" = names(components(subgraph_edges)$membership),"Subgraph" = components(subgraph_edges)$membership)

subgraph_membership %<>%
  arrange(Subgraph)

write_xlsx(subgraph_membership,"subgraph_membership_quant_20.xlsx")

# subgraph 1 contains the clade at the bottom of the tree
# this is ortho- and para-rubulaviruses including mumps, infecting primarily bats and primates (humans)

# subgraph 2 contains the cluster above the birds
# this comprises reptilian and fish-infecting viruses and a number of mammalian respiroviruses (mostly human and livestock)

# subgraph 3 contains the primarily carnivora + artiodactyla cluster
# this comprises mammalian-infecting morbilliviruses

# subgraph 4 contains the large top cluster
# this comprises mammalian-infecting jeilongviruses and henipaviruses

# subgraph 5 contains the cluster above the birds
# this comprises reptile-infecting viruses

# subgraph 6 contains a cluster between respiroviruses and morbilliviruses
# this comprises mammalian-infecting narmoviruses (almost all rodent viruses + one bovine)

# subgraph 7 is the orphan clade containing only 2 viruses, the Paraavulavirus wisconsinense
# these are avian infecting viruses on a very long branch from all the other avian viruses

# subgraph 8 contains the lower part of the avian-infecting clade
# this comprises the avian metaavulaviruses

# subgraph 9 contains the upper part of the avian-infecting clade
# this comprises the avian orthoavulaviruses

### manual subgraphs - find a threshold that preserves separate henipacluster with combined metaavulavirus cluster
quant <- quantile(igraph::E(g)$weight,probs = seq(0,1,0.01))

subgraph_edges <- c(igraph::E(g)[igraph::E(g)$weight<=33]) 
subgraph_edges <- subgraph_from_edges(g,subgraph_edges)

png(filename="./graph_phylo_plots/host_labels_subgraphs_final.png",width=1200,heigh=2400,units="px")
plot(subgraph_edges, vertex.size=3, vertex.label=NA,vertex.color=mycol)
dev.off()

# look at what sequences belong to which subclusters
subgraph_membership <- tibble("Accession" = names(components(subgraph_edges)$membership),"Subgraph" = components(subgraph_edges)$membership)

subgraph_membership %<>%
  arrange(Subgraph)

write_xlsx(subgraph_membership,"subgraph_membership_final.xlsx")
  
