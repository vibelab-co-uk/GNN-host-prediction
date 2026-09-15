library(tidyverse)
library(readxl)
library(magrittr)
library(plyr)
# this script generates simple feature sets using pre-generated features

gc <- read_csv("../../../4_feature_generation/3_genomic_feats/nucleotide_feats/paramyxo_gc_prop_features.csv")
nt <- read_csv("../../../4_feature_generation/3_genomic_feats/nucleotide_feats/paramyxo_nt_prop_features.csv")
dint <- read_csv("../../../4_feature_generation/3_genomic_feats/nucleotide_feats/paramyxo_dint_prop_features.csv")
kmer3 <- read_csv("../../../4_feature_generation/3_genomic_feats/nucleotide_feats/paramyxo_3mer_prop_features.csv")
kmer4 <- read_csv("../../../4_feature_generation/3_genomic_feats/nucleotide_feats/paramyxo_4mer_prop_features.csv")
kmer5 <- read_csv("../../../4_feature_generation/3_genomic_feats/nucleotide_feats/paramyxo_5mer_prop_features.csv")
kmer6 <- read_csv("../../../4_feature_generation/3_genomic_feats/nucleotide_feats/paramyxo_6mer_prop_features.csv")

rand1 <- read_csv("../../../4_feature_generation/3_genomic_feats/nucleotide_feats/paramyxo_rand_features_set1.csv")
rand2 <- read_csv("../../../4_feature_generation/3_genomic_feats/nucleotide_feats/paramyxo_rand_features_set2.csv")
rand3 <- read_csv("../../../4_feature_generation/3_genomic_feats/nucleotide_feats/paramyxo_rand_features_set3.csv")
rand4 <- read_csv("../../../4_feature_generation/3_genomic_feats/nucleotide_feats/paramyxo_rand_features_set4.csv")
rand5 <- read_csv("../../../4_feature_generation/3_genomic_feats/nucleotide_feats/paramyxo_rand_features_set5.csv")

aac <- read_csv("../../../4_feature_generation/3_genomic_feats/protein_feats/paramyxo_aac_prop_features.csv")
dp <- read_csv("../../../4_feature_generation/3_genomic_feats/protein_feats/paramyxo_dp_prop_features.csv")

new_gc <- read_csv("../4_feature_generation/3_genomic_feats/nucleotide_feats/paramyxo_gc_prop_features.csv")
new_nt <- read_csv("../4_feature_generation/3_genomic_feats/nucleotide_feats/paramyxo_nt_prop_features.csv")
new_dint <- read_csv("../4_feature_generation/3_genomic_feats/nucleotide_feats/paramyxo_dint_prop_features.csv")
new_kmer3 <- read_csv("../4_feature_generation/3_genomic_feats/nucleotide_feats/paramyxo_3mer_prop_features.csv")
new_kmer4 <- read_csv("../4_feature_generation/3_genomic_feats/nucleotide_feats/paramyxo_4mer_prop_features.csv")
new_kmer5 <- read_csv("../4_feature_generation/3_genomic_feats/nucleotide_feats/paramyxo_5mer_prop_features.csv")
new_kmer6 <- read_csv("../4_feature_generation/3_genomic_feats/nucleotide_feats/paramyxo_6mer_prop_features.csv")

new_rand1 <- read_csv("../4_feature_generation/3_genomic_feats/nucleotide_feats/paramyxo_rand_features_set1.csv")
new_rand2 <- read_csv("../4_feature_generation/3_genomic_feats/nucleotide_feats/paramyxo_rand_features_set2.csv")
new_rand3 <- read_csv("../4_feature_generation/3_genomic_feats/nucleotide_feats/paramyxo_rand_features_set3.csv")
new_rand4 <- read_csv("../4_feature_generation/3_genomic_feats/nucleotide_feats/paramyxo_rand_features_set4.csv")
new_rand5 <- read_csv("../4_feature_generation/3_genomic_feats/nucleotide_feats/paramyxo_rand_features_set5.csv")

new_aac <- read_csv("../4_feature_generation/3_genomic_feats/protein_feats/paramyxo_aac_prop_features.csv")
new_dp <- read_csv("../4_feature_generation/3_genomic_feats/protein_feats/paramyxo_dp_prop_features.csv")

full_gc <- rbind(gc,new_gc)
full_nt <- rbind(nt,new_nt)
full_dint <- rbind(dint,new_dint)
full_kmer3 <- rbind(kmer3,new_kmer3)
full_kmer4 <- rbind(kmer4,new_kmer4)
full_kmer5 <- rbind(kmer5,new_kmer5)
full_kmer6 <- rbind(kmer6,new_kmer6)

full_rand1 <- rbind(rand1,new_rand1)
full_rand2 <- rbind(rand2,new_rand2)
full_rand3 <- rbind(rand3,new_rand3)
full_rand4 <- rbind(rand4,new_rand4)
full_rand5 <- rbind(rand5,new_rand5)

full_aac <- rbind(aac,new_aac)
full_dp <- rbind(dp,new_dp)

master_list <- list(full_gc,full_nt,full_dint,full_aac,full_dp,full_kmer3,full_kmer4,full_kmer5,full_kmer6,full_rand1,full_rand2,full_rand3,full_rand4,full_rand5)

names(master_list) <- c("gc","nt","dint","aac","dp","kmer3","kmer4","kmer5","kmer6","rand1","rand2","rand3","rand4","rand5")

# read in table with cluster hosts
cluster_hosts <- read_csv("../3_host_assignment/combined_cluster_hosts.csv")
colnames(cluster_hosts)[1] <- "Accession"
cluster_hosts$Host_rank <- as.factor(cluster_hosts$Host_rank)

# basic feature sets
gc_only <- join_all(list(master_list["gc"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(gc_only,paste0("../6_GNN/feature_sets/count_feats_gconly.csv"))

dint_only <- join_all(list(master_list["dint"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(dint_only,paste0("../6_GNN/feature_sets/count_feats_dintonly.csv"))

aac_only <- join_all(list(master_list["aac"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(aac_only,paste0("../6_GNN/feature_sets/count_feats_aaconly.csv"))

dp_only <- join_all(list(master_list["dp"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(dp_only,paste0("../6_GNN/feature_sets/count_feats_dponly.csv"))

kmer3_only <- join_all(list(master_list["kmer3"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(kmer3_only,paste0("../6_GNN/feature_sets/count_feats_kmer3only.csv"))

# random pseudofeature sets
rand1 <- join_all(list(master_list["rand1"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(rand1,paste0("../6_GNN/feature_sets/count_feats_rand1.csv"))

rand2 <- join_all(list(master_list["rand2"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(rand2,paste0("../6_GNN/feature_sets/count_feats_rand2.csv"))

rand3 <- join_all(list(master_list["rand3"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(rand3,paste0("../6_GNN/feature_sets/count_feats_rand3.csv"))

rand4 <- join_all(list(master_list["rand4"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(rand4,paste0("../6_GNN/feature_sets/count_feats_rand4.csv"))

rand5 <- join_all(list(master_list["rand5"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(rand5,paste0("../6_GNN/feature_sets/count_feats_rand5.csv"))

