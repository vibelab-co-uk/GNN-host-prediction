library(tidyverse)
library(readxl)
library(magrittr)
library(plyr)

# this script generates simple feature sets using pre-generated features

setwd("/Users/jamieherzig/Documents/clean_run/4_feature_generation/3_genomic_feats/nucleotide_feats")

gc <- read_csv("paramyxo_gc_prop_features.csv")
nt <- read_csv("paramyxo_nt_prop_features.csv")
dint <- read_csv("paramyxo_dint_prop_features.csv")
kmer3 <- read_csv("paramyxo_3mer_prop_features.csv")
kmer4 <- read_csv("paramyxo_4mer_prop_features.csv")
kmer5 <- read_csv("paramyxo_5mer_prop_features.csv")
kmer6 <- read_csv("paramyxo_6mer_prop_features.csv")

setwd("/Users/jamieherzig/Documents/clean_run/4_feature_generation/3_genomic_feats/protein_feats")

aac <- read_csv("paramyxo_aac_prop_features.csv")
dp <- read_csv("paramyxo_dp_prop_features.csv")

master_list <- list(gc,nt,dint,aac,dp,kmer3,kmer4,kmer5,kmer6)

names(master_list) <- c("gc","nt","dint","aac","dp","kmer3","kmer4","kmer5","kmer6")

# read in table with cluster hosts
cluster_hosts <- read_csv("/Users/jamieherzig/Documents/clean_run/3_host_assignment/all_cluster_hosts.csv")
colnames(cluster_hosts)[1] <- "Accession"
cluster_hosts$Host_rank <- as.factor(cluster_hosts$Host_rank)

setwd("/Users/jamieherzig/Documents/clean_run/6_GNN/protein_ablation/simple_feats/feature_sets")

# basic feature sets
gc_only <- join_all(list(master_list["gc"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(gc_only,paste0("count_feats_gconly.csv"))

dint_only <- join_all(list(master_list["dint"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(dint_only,paste0("count_feats_dintonly.csv"))

aac_only <- join_all(list(master_list["aac"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(aac_only,paste0("count_feats_aaconly.csv"))

dp_only <- join_all(list(master_list["dp"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(dp_only,paste0("count_feats_dponly.csv"))

gc_aac <- join_all(list(master_list["gc"][[1]],master_list["aac"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(gc_aac,paste0("count_feats_gc_aac.csv"))

dint_dp <- join_all(list(master_list["dint"][[1]],master_list["dp"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(dint_dp,paste0("count_feats_dint_dp.csv"))

kmer3_only <- join_all(list(master_list["kmer3"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(kmer3_only,paste0("count_feats_kmer3only.csv"))

kmer3_dp <- join_all(list(master_list["kmer3"][[1]],master_list["dp"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(kmer3_dp,paste0("count_feats_kmer3_dp.csv"))

all_basic <- join_all(list(master_list["gc"][[1]],master_list["nt"][[1]],master_list["dint"][[1]],master_list["aac"][[1]],master_list["dp"][[1]],
                           cluster_hosts), by='Accession', type='inner')
write_csv(all_basic,paste0("count_feats_allbasic.csv"))

# generate a dummy kmer feature set where every kmer is equally likely
dummy_kmer <- kmer3_only

dummy_kmer[,2:(ncol(dummy_kmer)-1)] <- 1/64 # 64 3mers
write_csv(dummy_kmer,paste0("count_feats_dummykmer.csv"))

# finally write a text file containing folder names to generate results folders
namelist <- vector()

for(i in 1:nrow(all_combs)){
  namelist[i] <- paste0("results_",all_combs[i,1],"_",all_combs[i,2])
}

namelist[length(namelist)+1] <- paste0("results_gconly")
namelist[length(namelist)+1] <- paste0("results_dintonly")
namelist[length(namelist)+1] <- paste0("results_aaconly")
namelist[length(namelist)+1] <- paste0("results_dponly")
namelist[length(namelist)+1] <- paste0("results_gc_aac")
namelist[length(namelist)+1] <- paste0("results_dint_dp")
namelist[length(namelist)+1] <- paste0("results_kmer3only")
namelist[length(namelist)+1] <- paste0("results_kmer3_dp")
namelist[length(namelist)+1] <- paste0("results_allbasic")
namelist[length(namelist)+1] <- paste0("results_dummykmer")

# print list of directory names 1 per line for bash scripting
cat(namelist,file="dir_names.txt",sep="
")
