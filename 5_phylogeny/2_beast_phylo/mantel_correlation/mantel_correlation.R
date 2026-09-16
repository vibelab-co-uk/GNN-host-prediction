
#Load packages------------------------------------------------------------------

library(ape)
library(vegan)
library(dplyr)
library(tidyr)
library(readr)


#Load data----------------------------------------------------------------------

tree_file <- "GTR_HMC_HPSTR_tree_50_pruned.tree"
host_order_matrix <- read_csv("all_cluster_hosts.csv")


#Read tree ---------------------------------------------------------------------

tree <- read.tree(tree_file)


#Create phylogenetic distance matrix--------------------------------------------

phylo_dist <- cophenetic(tree)


#Create binary host rank matrix-------------------------------------------------

host_order_matrix <- host_order_matrix |>
  mutate(value = 1) |>
  pivot_wider(
    names_from = Host_rank,
    values_from = value,
    values_fill = 0)


#Accession column to rowname---------------------------------------------------

host_order_matrix <- as.data.frame(host_order_matrix)

rownames(host_order_matrix) <- host_order_matrix$ref_accessions

host_order_matrix <- host_order_matrix |> 
  select(-ref_accessions)


View(host_order_matrix)


#Check that all tree tips are present-------------------------------------------

stopifnot(all(rownames(phylo_dist) %in% rownames(host_order_matrix)))


#Reorder rows to match phylogenetic matrix-------------------------------------

host_order_matrix <- host_order_matrix[
  rownames(phylo_dist),
  ,
  drop = FALSE]


#Create trait distance matrix---------------------------------------------------

host_dist_order <- vegdist(
  host_order_matrix,
  method = "jaccard",
  binary = T)


#Check matrices are aligned-----------------------------------------------------

identical(
  rownames(phylo_dist),
  attr(host_dist_order, "Labels"))


#Run Mantel test----------------------------------------------------------------

set.seed(57)

mantel_result <- mantel(
  phylo_dist,
  host_dist_order,
  method = "pearson",
  permutations = 9999)


mantel_result


