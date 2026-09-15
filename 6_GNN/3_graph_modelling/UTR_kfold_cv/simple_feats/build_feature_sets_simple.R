library(tidyverse)
library(readxl)
library(magrittr)
library(plyr)

# this script generates simple feature sets using pre-generated features
# both leader & trailer UTRs ----------------------------------------------
gc_both <- read_csv("../../../../4_feature_generation/3_genomic_feats/UTR_feats/paramyxo_gc_prop_bothUTR_features.csv")
nt_both <- read_csv("../../../../4_feature_generation/3_genomic_feats/UTR_feats/paramyxo_nt_prop_bothUTR_features.csv")
dint_both <- read_csv("../../../../4_feature_generation/3_genomic_feats/UTR_feats/paramyxo_dint_prop_bothUTR_features.csv")
kmer3_both <- read_csv("../../../../4_feature_generation/3_genomic_feats/UTR_feats/paramyxo_3mer_prop_bothUTR_features.csv")
kmer4_both <- read_csv("../../../../4_feature_generation/3_genomic_feats/UTR_feats/paramyxo_4mer_prop_bothUTR_features.csv")
kmer5_both <- read_csv("../../../../4_feature_generation/3_genomic_feats/UTR_feats/paramyxo_5mer_prop_bothUTR_features.csv")
kmer6_both <- read_csv("../../../../4_feature_generation/3_genomic_feats/UTR_feats/paramyxo_6mer_prop_bothUTR_features.csv")

master_list <- list(gc_both,nt_both,dint_both,kmer3_both,kmer4_both,kmer5_both,kmer6_both)

names(master_list) <- c("gc_both","nt_both","dint_both","kmer3_both","kmer4_both","kmer5_both","kmer6_both")

# read in table with cluster hosts
cluster_hosts <- read_csv("../../../../3_host_assignment/all_cluster_hosts.csv")
colnames(cluster_hosts)[1] <- "Accession"
cluster_hosts$Host_rank <- as.factor(cluster_hosts$Host_rank)

# basic feature sets
gc_only <- join_all(list(master_list["gc_both"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(gc_only,paste0("./feature_sets/count_feats_both_gconly.csv"))

dint_only <- join_all(list(master_list["dint_both"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(dint_only,paste0("./feature_sets/count_feats_both_dintonly.csv"))

kmer3_only <- join_all(list(master_list["kmer3_both"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(kmer3_only,paste0("./feature_sets/count_feats_both_kmer3only.csv"))

kmer4_only <- join_all(list(master_list["kmer4_both"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(kmer4_only,paste0("./feature_sets/count_feats_both_kmer4only.csv"))

kmer5_only <- join_all(list(master_list["kmer5_both"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(kmer5_only,paste0("./feature_sets/count_feats_both_kmer5only.csv"))

kmer6_only <- join_all(list(master_list["kmer6_both"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(kmer6_only,paste0("./feature_sets/count_feats_both_kmer6only.csv"))


# leader_only ----------------------------------------------
gc_leader <- read_csv("../../../../4_feature_generation/3_genomic_feats/UTR_feats/paramyxo_gc_prop_leaderUTR_features.csv")
nt_leader <- read_csv("../../../../4_feature_generation/3_genomic_feats/UTR_feats/paramyxo_nt_prop_leaderUTR_features.csv")
dint_leader <- read_csv("../../../../4_feature_generation/3_genomic_feats/UTR_feats/paramyxo_dint_prop_leaderUTR_features.csv")
kmer3_leader <- read_csv("../../../../4_feature_generation/3_genomic_feats/UTR_feats/paramyxo_3mer_prop_leaderUTR_features.csv")
kmer4_leader <- read_csv("../../../../4_feature_generation/3_genomic_feats/UTR_feats/paramyxo_4mer_prop_leaderUTR_features.csv")
kmer5_leader <- read_csv("../../../../4_feature_generation/3_genomic_feats/UTR_feats/paramyxo_5mer_prop_leaderUTR_features.csv")
kmer6_leader <- read_csv("../../../../4_feature_generation/3_genomic_feats/UTR_feats/paramyxo_6mer_prop_leaderUTR_features.csv")

master_list <- list(gc_leader,nt_leader,dint_leader,kmer3_leader,kmer4_leader,kmer5_leader,kmer6_leader)

names(master_list) <- c("gc_leader","nt_leader","dint_leader","kmer3_leader","kmer4_leader","kmer5_leader","kmer6_leader")

# read in table with cluster hosts
cluster_hosts <- read_csv("../../../../3_host_assignment/all_cluster_hosts.csv")
colnames(cluster_hosts)[1] <- "Accession"
cluster_hosts$Host_rank <- as.factor(cluster_hosts$Host_rank)

# basic feature sets
gc_only <- join_all(list(master_list["gc_leader"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(gc_only,paste0("./feature_sets/count_feats_leader_gconly.csv"))

dint_only <- join_all(list(master_list["dint_leader"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(dint_only,paste0("./feature_sets/count_feats_leader_dintonly.csv"))

kmer3_only <- join_all(list(master_list["kmer3_leader"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(kmer3_only,paste0("./feature_sets/count_feats_leader_kmer3only.csv"))

kmer4_only <- join_all(list(master_list["kmer4_leader"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(kmer4_only,paste0("./feature_sets/count_feats_leader_kmer4only.csv"))

kmer5_only <- join_all(list(master_list["kmer5_leader"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(kmer5_only,paste0("./feature_sets/count_feats_leader_kmer5only.csv"))

kmer6_only <- join_all(list(master_list["kmer6_leader"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(kmer6_only,paste0("./feature_sets/count_feats_leader_kmer6only.csv"))

# trailer only ----------------------------------------------
gc_trailer <- read_csv("../../../../4_feature_generation/3_genomic_feats/UTR_feats/paramyxo_gc_prop_trailerUTR_features.csv")
nt_trailer <- read_csv("../../../../4_feature_generation/3_genomic_feats/UTR_feats/paramyxo_nt_prop_trailerUTR_features.csv")
dint_trailer <- read_csv("../../../../4_feature_generation/3_genomic_feats/UTR_feats/paramyxo_dint_prop_trailerUTR_features.csv")
kmer3_trailer <- read_csv("../../../../4_feature_generation/3_genomic_feats/UTR_feats/paramyxo_3mer_prop_trailerUTR_features.csv")
kmer4_trailer <- read_csv("../../../../4_feature_generation/3_genomic_feats/UTR_feats/paramyxo_4mer_prop_trailerUTR_features.csv")
kmer5_trailer <- read_csv("../../../../4_feature_generation/3_genomic_feats/UTR_feats/paramyxo_5mer_prop_trailerUTR_features.csv")
kmer6_trailer <- read_csv("../../../../4_feature_generation/3_genomic_feats/UTR_feats/paramyxo_6mer_prop_trailerUTR_features.csv")

master_list <- list(gc_trailer,nt_trailer,dint_trailer,kmer3_trailer,kmer4_trailer,kmer5_trailer,kmer6_trailer)

names(master_list) <- c("gc_trailer","nt_trailer","dint_trailer","kmer3_trailer","kmer4_trailer","kmer5_trailer","kmer6_trailer")

# read in table with cluster hosts
cluster_hosts <- read_csv("../../../../3_host_assignment/all_cluster_hosts.csv")
colnames(cluster_hosts)[1] <- "Accession"
cluster_hosts$Host_rank <- as.factor(cluster_hosts$Host_rank)

# basic feature sets
gc_only <- join_all(list(master_list["gc_trailer"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(gc_only,paste0("./feature_sets/count_feats_trailer_gconly.csv"))

dint_only <- join_all(list(master_list["dint_trailer"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(dint_only,paste0("./feature_sets/count_feats_trailer_dintonly.csv"))

kmer3_only <- join_all(list(master_list["kmer3_trailer"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(kmer3_only,paste0("./feature_sets/count_feats_trailer_kmer3only.csv"))

kmer4_only <- join_all(list(master_list["kmer4_trailer"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(kmer4_only,paste0("./feature_sets/count_feats_trailer_kmer4only.csv"))

kmer5_only <- join_all(list(master_list["kmer5_trailer"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(kmer5_only,paste0("./feature_sets/count_feats_trailer_kmer5only.csv"))

kmer6_only <- join_all(list(master_list["kmer6_trailer"][[1]],cluster_hosts), by='Accession', type='inner')
write_csv(kmer6_only,paste0("./feature_sets/count_feats_trailer_kmer6only.csv"))

# finally write a text file containing folder names to generate results folders
namelist <- vector()

namelist[length(namelist)+1] <- paste0("results_both_gconly")
namelist[length(namelist)+1] <- paste0("results_both_dintonly")
namelist[length(namelist)+1] <- paste0("results_both_kmer3only")
namelist[length(namelist)+1] <- paste0("results_both_kmer4only")
namelist[length(namelist)+1] <- paste0("results_both_kmer5only")
namelist[length(namelist)+1] <- paste0("results_both_kmer6only")
namelist[length(namelist)+1] <- paste0("results_leader_gconly")
namelist[length(namelist)+1] <- paste0("results_leader_dintonly")
namelist[length(namelist)+1] <- paste0("results_leader_kmer3only")
namelist[length(namelist)+1] <- paste0("results_leader_kmer4only")
namelist[length(namelist)+1] <- paste0("results_leader_kmer5only")
namelist[length(namelist)+1] <- paste0("results_leader_kmer6only")
namelist[length(namelist)+1] <- paste0("results_trailer_gconly")
namelist[length(namelist)+1] <- paste0("results_trailer_dintonly")
namelist[length(namelist)+1] <- paste0("results_trailer_kmer3only")
namelist[length(namelist)+1] <- paste0("results_trailer_kmer4only")
namelist[length(namelist)+1] <- paste0("results_trailer_kmer5only")
namelist[length(namelist)+1] <- paste0("results_trailer_kmer6only")
# print list of directory names 1 per line for bash scripting
cat(namelist,file="./feature_sets/dir_names.txt",sep="
")
