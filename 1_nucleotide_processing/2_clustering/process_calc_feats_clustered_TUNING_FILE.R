library(tidyverse)
library(rentrez)
library(Biostrings)
library(taxizedb)
library(coRdon)
library(protr)
library(magrittr)
library(writexl)

# read and process nucleotide metadata-----------------------------------------------------
tune_values <- read_tsv("current_values.txt",col_names = F)
setwd(paste0(tune_values[1,1],"_",tune_values[2,1],"_clusters"))

# read in metadata table including the host information for all paramyxoviruses
paramyxo_nuc_metadata <- read_csv("../1_initial_cleaning/11158_paramyxoviridae_complete_metadata_21112025.csv")

# drop entries without host data
paramyxo_nuc_metadata <- paramyxo_nuc_metadata %>%
  drop_na(Host)

# now we want to find the parent taxonomy at given level for every host
nuc_host_taxa <- name2taxid(paramyxo_nuc_metadata$Host,db="ncbi",out_type="summary") # convert taxon names to IDs
nuc_host_taxa$rank_Class <- NA
nuc_host_taxa$rank_Order <- NA

nuc_full_class_list <- classification(nuc_host_taxa$id)
# 
# taxassign <- function(tax_id){
#   rsp <- deparse(substitute(tax_id))
#   
#   nuc_host_taxa <- nuc_host_taxa %>% 
#     mutate(rank_Class = case_when(id == !!tax_id ~ (full_class_list[[rsp]] %>% filter(rank=="class")) %>% pull(name),
#            .default=rank_Class))
#   return(nuc_host_taxa)
# }
# 
# lapply(as.numeric(nuc_host_taxa$id),taxassign)

for(i in 1:nrow(nuc_host_taxa)){
  if(length(nuc_full_class_list[i][[1]] %>% filter(rank=="class") %>% pull(name))==0){
  nuc_host_taxa[i,3] <- nuc_full_class_list[i][[1]] %>% filter(rank=="superclass") %>% pull(name)
  } else {
    nuc_host_taxa[i,3] <- nuc_full_class_list[i][[1]] %>% filter(rank=="class") %>% pull(name)
  }
}

for(i in 1:nrow(nuc_host_taxa)){
  if((length(nuc_full_class_list[i][[1]] %>% filter(rank=="superorder") %>% pull(name))==0) &
      (length(nuc_full_class_list[i][[1]] %>% filter(rank=="order") %>% pull(name))==0)) {
    nuc_host_taxa[i,4] <- NA
  } else if(length(nuc_full_class_list[i][[1]] %>% filter(rank=="order") %>% pull(name))==0){
    nuc_host_taxa[i,4] <- nuc_full_class_list[i][[1]] %>% filter(rank=="superorder") %>% pull(name)
  } else {
    nuc_host_taxa[i,4] <- nuc_full_class_list[i][[1]] %>% filter(rank=="order") %>% pull(name)
  }
}

# check if any of the taxa result in a duplication

if(length(unique(nuc_host_taxa$name))!=nrow(nuc_host_taxa)){
  print(nuc_host_taxa$name[duplicated(nuc_host_taxa$name)]
)
}

# manually correct these
nuc_host_taxa <- nuc_host_taxa %>% 
  filter(rank_Class != "Magnoliopsida")

# join this back to the metadata
colnames(nuc_host_taxa) <- c("Host","Host_id","Host_rank_class","Host_rank_order")

paramyxo_nuc_metadata <- inner_join(paramyxo_nuc_metadata,nuc_host_taxa,by="Host")

# read and process nucleotide sequences ----------------------------------------------
# read in a TSV containing list of sequences for each host
clusters <- read.delim("output.FASTA_cluster.tsv",header=F)

# for each cluster, get the accessions for every member
colnames(clusters) <- c("ref_accessions","members")

# assign host taxID using the cleaned IDs in the metadata
host_IDs <- paramyxo_nuc_metadata %>% 
  dplyr::select(Host_id,Accession)

clusters <- left_join(clusters,host_IDs,join_by(ref_accessions == Accession))
clusters <- left_join(clusters,host_IDs,join_by(members == Accession))

organism_name <- paramyxo_nuc_metadata %>% 
  dplyr::select(Organism_Name,Accession)

clusters <- left_join(clusters,organism_name,join_by(ref_accessions == Accession))
clusters <- left_join(clusters,organism_name,join_by(members == Accession))

colnames(clusters) <- c("ref_accessions","members","ref_host","member_host","ref_organism","member_organism")

# use stringr to try and do some basic cleaning of virus species names
clusters$ref_organism <- stringr::str_to_sentence(clusters$ref_organism)
clusters$member_organism <- stringr::str_to_sentence(clusters$member_organism)

# do some manual cleaning of organism names to deal with 
clusters %<>%
  mutate(ref_organism = case_when(ref_organism == "Paraavulavirus hongkongense" ~ "Avian paramyxovirus 4",
                                  ref_organism == "Metaavulavirus hongkongense" ~ "Avian metaavulavirus 6",
                                  ref_organism == "Avian paramyxovirus 6" ~ "Avian metaavulavirus 6",
                                  ref_organism == "Avian paramyxovirus penguin/falkland islands/324/2007" ~ "Avian paramyxovirus 10",
                                  ref_organism == "Jeilongvirus beilongi" ~ "Beilong virus",
                                  grepl("Rinderpest",ref_organism,ignore.case = T) ~ "Rinderpest morbillivirus",
                                  ref_organism == "Canine distemper virus" ~ "Morbillivirus canis",
                                  ref_organism == "Human parainfluenza virus 4a" ~ "Human orthorubulavirus 4",
                                  ref_organism == "Human parainfluenza virus 4b" ~ "Human orthorubulavirus 4",
                                  ref_organism == "Small ruminant morbillivirus" ~ "Morbillivirus caprinae",
                                  grepl("Peste",ref_organism,ignore.case = T) ~ "Morbillivirus caprinae",
                                  ref_organism == "Avian paramyxovirus 1" ~ "Avian orthoavulavirus 1",
                                  ref_organism == "Orthoavulavirus javaense" ~ "Avian orthoavulavirus 1",
                                  ref_organism == "Pigeon paramyxovirus 1" ~ "Avian orthoavulavirus 1",
                                  grepl("Newcastle|NDV",ref_organism,ignore.case = T) ~ "Avian orthoavulavirus 1",
                                  grepl("Measles",ref_organism,ignore.case = T) ~ "Measles morbillivirus",
                                  ref_organism == "Morbillivirus hominis" ~ "Measles morbillivirus",
                                  grepl("Feline morbillivirus",ref_organism,ignore.case = T) ~ "Feline morbillivirus",
                                  ref_organism == "Orthorubulavirus laryngotracheitidis" ~ "Human orthorubulavirus 2",
                                  grepl("Avian paramyxovirus 13",ref_organism,ignore.case = T) ~ "Avian orthoavulavirus 13",
                                  ref_organism == "Respirovirus laryngotracheitidis" ~ "Human respirovirus 1",
                                  ref_organism == "Dolphin morbillivirus" ~ "Morbillivirus ceti",
                                  ref_organism == "Ovine parainfluenza virus 3" ~ "Caprine parainfluenza virus 3",
                                  ref_organism == "Respirovirus pneumoniae" ~ "Human respirovirus 3",
                                  ref_organism == "Parainfluenza virus 5" ~ "Mammalian orthorubulavirus 5",
                                  ref_organism == "Orthoavulavirus newyorkense" ~ "Avian orthoavulavirus 9",
                                  ref_organism == "Avian paramyxovirus 10" ~ "Avian metaavulavirus 10",
                                  ref_organism == "Metaavulavirus peixense" ~ "Avian paramyxovirus 15",
                                  grepl("Mumps",ref_organism,ignore.case = T) ~ "Mumps orthorubulavirus",
                                  .default = ref_organism)) %>% 
  mutate(member_organism = case_when(member_organism == "Paraavulavirus hongkongense" ~ "Avian paramyxovirus 4",
                                  member_organism == "Metaavulavirus hongkongense" ~ "Avian metaavulavirus 6",
                                  member_organism == "Avian paramyxovirus 6" ~ "Avian metaavulavirus 6",
                                  member_organism == "Avian paramyxovirus penguin/falkland islands/324/2007" ~ "Avian paramyxovirus 10",
                                  member_organism == "Jeilongvirus beilongi" ~ "Beilong virus",
                                  grepl("Rinderpest",member_organism,ignore.case = T) ~ "Rinderpest morbillivirus",
                                  member_organism == "Canine distemper virus" ~ "Morbillivirus canis",
                                  member_organism == "Human parainfluenza virus 4a" ~ "Human orthorubulavirus 4",
                                  member_organism == "Human parainfluenza virus 4b" ~ "Human orthorubulavirus 4",
                                  member_organism == "Small ruminant morbillivirus" ~ "Morbillivirus caprinae",
                                  grepl("Peste",member_organism,ignore.case = T) ~ "Morbillivirus caprinae",
                                  member_organism == "Avian paramyxovirus 1" ~ "Avian orthoavulavirus 1",
                                  member_organism == "Orthoavulavirus javaense" ~ "Avian orthoavulavirus 1",
                                  member_organism == "Pigeon paramyxovirus 1" ~ "Avian orthoavulavirus 1",
                                  grepl("Newcastle|NDV",member_organism,ignore.case = T) ~ "Avian orthoavulavirus 1",
                                  grepl("Measles",member_organism,ignore.case = T) ~ "Measles morbillivirus",
                                  member_organism == "Morbillivirus hominis" ~ "Measles morbillivirus",
                                  grepl("Feline morbillivirus",member_organism,ignore.case = T) ~ "Feline morbillivirus",
                                  member_organism == "Orthorubulavirus laryngotracheitidis" ~ "Human orthorubulavirus 2",
                                  grepl("Avian paramyxovirus 13",member_organism,ignore.case = T) ~ "Avian orthoavulavirus 13",
                                  member_organism == "Respirovirus laryngotracheitidis" ~ "Human respirovirus 1",
                                  member_organism == "Dolphin morbillivirus" ~ "Morbillivirus ceti",
                                  member_organism == "Ovine parainfluenza virus 3" ~ "Caprine parainfluenza virus 3",
                                  member_organism == "Respirovirus pneumoniae" ~ "Human respirovirus 3",
                                  member_organism == "Parainfluenza virus 5" ~ "Mammalian orthorubulavirus 5",
                                  member_organism == "Orthoavulavirus newyorkense" ~ "Avian orthoavulavirus 9",
                                  member_organism == "Avian paramyxovirus 10" ~ "Avian metaavulavirus 10",
                                  member_organism == "Metaavulavirus peixense" ~ "Avian paramyxovirus 15",
                                  grepl("Mumps",member_organism,ignore.case = T) ~ "Mumps orthorubulavirus",
                                  .default = member_organism))

# generate a list of all recorded hosts for every cluster
cluster_host_search <- function(x){
  
  indiv_cluster <- clusters %>% 
    filter(ref_accessions == x) 
  
  indiv_cluster <- tibble(ref_accession = NA, hosts = unique(c(indiv_cluster$ref_host,indiv_cluster$member_host)))
  
  indiv_cluster$ref_accession <- x
  
  indiv_cluster %<>%
    drop_na()
  
  return(indiv_cluster)
}

all_cluster_hosts <- do.call(rbind,lapply(unique(clusters$ref_accessions),cluster_host_search))

tax_info <- paramyxo_nuc_metadata %>% 
  dplyr::select(Host_id,Host_rank_order,Host_rank_class) %>% 
  distinct()

all_cluster_hosts <- left_join(all_cluster_hosts,tax_info,join_by(hosts == Host_id))

# check if any clusters have no host information - if so, drop them
nrow(all_cluster_hosts %>% 
       filter(is.na(hosts)==T))

all_cluster_hosts %<>% 
  drop_na(hosts)

order_raw <- all_cluster_hosts %>% 
  dplyr::count(ref_accession,Host_rank_order)

# this double count first sums the number of taxonomic orders for each ref accession and then counts the numbers of entries for each accession
# therefore count > 1 implies the cluster has members with hosts from different taxonomic orders
order_test <- all_cluster_hosts %>% 
  dplyr::count(ref_accession,Host_rank_order) %>% 
  dplyr::count(ref_accession) %>% 
  filter(n>1)

nrow(order_test)

class_raw <- all_cluster_hosts %>% 
  dplyr::count(ref_accession,Host_rank_class)

# this double count first sums the number of taxonomic classes for each ref accession and then counts the numbers of entries for each accession
# therefore count > 1 implies the cluster has members with hosts from different taxonomic classes
class_test <- all_cluster_hosts %>% 
  dplyr::count(ref_accession,Host_rank_class) %>% 
  dplyr::count(ref_accession) %>% 
  filter(n>1)

nrow(class_test)

# generate a list of all recorded viruses for every cluster
cluster_organism_search <- function(x){
  
  indiv_cluster <- clusters %>% 
    filter(ref_accessions == x) 
  
  indiv_cluster <- tibble(ref_accession = NA, organism = c(indiv_cluster$ref_organism,indiv_cluster$member_organism))
  
  indiv_cluster$ref_accession <- x
  
  indiv_cluster %<>%
    drop_na()
  
  return(indiv_cluster)
}


all_cluster_organism <- do.call(rbind,lapply(unique(clusters$ref_accessions),cluster_organism_search))

organism_test <- all_cluster_organism %>% 
  dplyr::count(ref_accession,organism)

# this double count first calculates the number of instances of each viral species in a cluster and then looks for clusters which only contain a single viral species
organism_summary <- all_cluster_organism %>% 
  dplyr::count(ref_accession,organism) %>% 
  dplyr::count(ref_accession) %>% 
  filter(n==1)

nrow(organism_summary)

# alternate count which instead records the number of clusters representing n viral species
# we'll use this to plot the distribution of clusters containing different numbers of viral species
multi_organism_test <- all_cluster_organism %>% 
  dplyr::count(ref_accession,organism) %>% 
  dplyr::count(ref_accession) %>% 
  dplyr::count(n)

colnames(multi_organism_test) <- c("species","clusters")

# calculate as a proportion of the number of clusters
multi_organism_test %<>%
  mutate(cluster_prop = clusters/(length(unique(all_cluster_hosts$ref_accession))))

# finally calculate how many cluster centroids represent the same species
centroid_test <- clusters %>% 
  drop_na() %>% 
  dplyr::distinct(ref_accessions,ref_organism)

centroid_test %<>%
  dplyr::count(ref_organism)

centroid_summary <- centroid_test %>%
  dplyr::count(n)

colnames(centroid_summary) <-c("number_of_centroids","count")

# calculate as a proportion of the number of clusters
centroid_summary %<>%
  mutate(count_prop = count/(length(unique(all_cluster_hosts$ref_accession))))

# calculate some basic statistics about cluster sizes and singleton clusters
cluster_sizes <- clusters %>%
  group_by(ref_accessions) %>%
  summarise(
    cluster_size = n(),
    .groups = 'drop'
  )

size_metrics <- cluster_sizes %>% 
  summarise(
    total_clusters = n(),
    singleton_clusters = sum(cluster_size == 1, na.rm = TRUE),
    percentage_singletons = (total_clusters - singleton_clusters) / total_clusters,
    mean_cluster_size = mean(cluster_size, na.rm = TRUE),
    median_cluster_size = median(cluster_size, na.rm = TRUE),
    .groups = 'drop'
  )

output <- tibble("n_clusters"=length(unique(clusters$ref_accessions)),
                 "multi_order_prop"=nrow(order_test)/length(unique(all_cluster_hosts$ref_accession)),
                 "multi_class_prop"=nrow(class_test)/length(unique(all_cluster_hosts$ref_accession)),
                 "mono_viral_species_prop"=length(unique(organism_test$ref_accession))/length(unique(all_cluster_organism$ref_accession)))

write_csv(output,paste0("min_seq_id_",tune_values[1,1],"_coverage_",tune_values[2,1],"_stats.csv"))
write_csv(multi_organism_test,paste0("min_seq_id_",tune_values[1,1],"_coverage_",tune_values[2,1],"_virus_species_cluster_count.csv"))
write_csv(organism_test,paste0("min_seq_id_",tune_values[1,1],"_coverage_",tune_values[2,1],"_virus_species_raw.csv"))
write_csv(order_raw,paste0("min_seq_id_",tune_values[1,1],"_coverage_",tune_values[2,1],"_orders.csv"))
write_csv(class_raw,paste0("min_seq_id_",tune_values[1,1],"_coverage_",tune_values[2,1],"_classes.csv"))
write_csv(centroid_test,paste0("min_seq_id_",tune_values[1,1],"_coverage_",tune_values[2,1],"_centroid_species_raw.csv"))
write_csv(centroid_summary,paste0("min_seq_id_",tune_values[1,1],"_coverage_",tune_values[2,1],"_centroid_species_cluster_count.csv"))
write_csv(size_metrics,paste0("min_seq_id_",tune_values[1,1],"_coverage_",tune_values[2,1],"_cluster_size_metrics.csv"))
write_csv(cluster_sizes,paste0("min_seq_id_",tune_values[1,1],"_coverage_",tune_values[2,1],"_cluster_size.csv"))
