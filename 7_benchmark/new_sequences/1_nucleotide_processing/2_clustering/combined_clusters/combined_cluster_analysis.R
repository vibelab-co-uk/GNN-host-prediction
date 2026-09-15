library(tidyverse)
library(Biostrings)

`%nin%` <- Negate(`%in%`)

combined_clusters <- read_tsv("clusters_09_04_combined.FASTA_cluster.tsv",col_names = c("Centroid","Members"))
new_clusters <- read_tsv("../clusters_09_04.FASTA_cluster.tsv",col_names = c("Centroid","Members"))
old_clusters <- read_tsv("../../../../../1_nucleotide_processing/2_clustering/0.9_0.4_clusters/output.FASTA_cluster.tsv",col_names = c("Centroid","Members"))

new_seq_ids <- unique(new_clusters$Centroid)

combined_clusters %<>% 
  mutate(new_seq = case_when(Members %in% (new_seq_ids) == T ~ T,
                             .default = F))

combined_clusters %<>% 
  group_by(Centroid) %>% 
  mutate(mix = case_when(length(unique(new_seq)) > 1 ~ T,
                         .default = F)) %>% 
  mutate(composition = case_when(new_seq %in% T & mix == F ~ "All new",
                                 new_seq %in% F & mix == F ~ "All old",
                                 .default = "Mixed cluster"))

cluster_counts <- combined_clusters %>% 
  ungroup() %>% 
  distinct(Centroid,composition) %>% 
  dplyr::count(composition)

# new sequences form 9 brand new clusters 
write_csv(cluster_counts, "cluster_counts.csv")
