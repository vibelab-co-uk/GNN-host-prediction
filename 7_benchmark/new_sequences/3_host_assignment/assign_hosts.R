library(tidyverse)
library(rentrez)
library(Biostrings)
library(taxizedb)
library(coRdon)
library(protr)
library(magrittr)
library(writexl)
library(readxl)
library(fishualize)

`%nin%` <- Negate(`%in%`)

# read and process nucleotide metadata -----------------------------------------------------
# read in list validated virus host lookup table - only viral species present in new dataset and not old need to be validated
lookup_tbl <- read_csv("../../../3_host_assignment/virus_host_lookup.csv")

# # read in clean metadata of prior sequences
# # this was used for confirming novel sequences in new data
# old_metadata <- read_xlsx("../../../3_host_assignment/paramyxovirus_metadata_cleaned_organism_names.xlsx")

# read in metadata table including the host information for all paramyxoviruses
paramyxo_nuc_metadata <- read_csv("../1_nucleotide_processing/1_initial_cleaning/1_quick_clean/11158_paramyxoviridae_complete_metadata_05112025_17062026.csv")

# clean virus species names
paramyxo_nuc_metadata$Organism_Name <- stringr::str_to_sentence(paramyxo_nuc_metadata$Organism_Name)

paramyxo_nuc_metadata %<>%
  mutate(Organism_Name = case_when(Organism_Name == "Orthoavulavirus upoense" ~ "Avian orthoavulavirus 16", # new correction required for new sequences
                                   Organism_Name == "Orthorubulavirus parotitidis" ~ "Mumps orthorubulavirus", # new correction required for new sequences
                                   Organism_Name == "Avian paramyxovirus 9" ~ "Avian orthoavulavirus 9", # new correction required for new sequences
                                   Organism_Name == "Paraavulavirus hongkongense" ~ "Avian paramyxovirus 4",
                                   Organism_Name == "Metaavulavirus hongkongense" ~ "Avian metaavulavirus 6",
                                   Organism_Name == "Avian paramyxovirus 6" ~ "Avian metaavulavirus 6",
                                   Organism_Name == "Avian paramyxovirus penguin/falkland islands/324/2007" ~ "Avian paramyxovirus 10",
                                   Organism_Name == "Jeilongvirus beilongi" ~ "Beilong virus",
                                   grepl("Rinderpest",Organism_Name,ignore.case = T) ~ "Rinderpest morbillivirus",
                                   grepl("Canine distemper virus",Organism_Name,ignore.case = T) ~ "Morbillivirus canis",
                                   Organism_Name == "Human parainfluenza virus 4a" ~ "Human orthorubulavirus 4",
                                   Organism_Name == "Human parainfluenza virus 4b" ~ "Human orthorubulavirus 4",
                                   grepl("Human parainfluenza virus 1",Organism_Name,ignore.case = T) ~ "Human respirovirus 1",
                                   Organism_Name == "Small ruminant morbillivirus" ~ "Morbillivirus caprinae",
                                   grepl("Peste",Organism_Name,ignore.case = T) ~ "Morbillivirus caprinae",
                                   Organism_Name == "Avian paramyxovirus 1" ~ "Avian orthoavulavirus 1",
                                   Organism_Name == "Orthoavulavirus javaense" ~ "Avian orthoavulavirus 1",
                                   Organism_Name == "Pigeon paramyxovirus 1" ~ "Avian orthoavulavirus 1",
                                   Organism_Name == "Goose paramyxovirus sf02" ~ "Avian orthoavulavirus 1",
                                   grepl("Newcastle|NDV",Organism_Name,ignore.case = T) ~ "Avian orthoavulavirus 1",
                                   grepl("Measles",Organism_Name,ignore.case = T) ~ "Measles morbillivirus",
                                   Organism_Name == "Morbillivirus hominis" ~ "Measles morbillivirus",
                                   grepl("Feline morbillivirus",Organism_Name,ignore.case = T) ~ "Feline morbillivirus",
                                   Organism_Name == "Orthorubulavirus laryngotracheitidis" ~ "Human orthorubulavirus 2",
                                   grepl("Avian paramyxovirus 13",Organism_Name,ignore.case = T) ~ "Avian orthoavulavirus 13",
                                   Organism_Name == "Respirovirus laryngotracheitidis" ~ "Human respirovirus 1",
                                   Organism_Name == "Dolphin morbillivirus" ~ "Morbillivirus ceti",
                                   Organism_Name == "Ovine parainfluenza virus 3" ~ "Caprine parainfluenza virus 3",
                                   Organism_Name == "Respirovirus pneumoniae" ~ "Human respirovirus 3",
                                   Organism_Name == "Parainfluenza virus 5" ~ "Mammalian orthorubulavirus 5",
                                   Organism_Name == "Orthorubulavirus mammalis" ~ "Mammalian orthorubulavirus 5",
                                   Organism_Name == "Orthoavulavirus newyorkense" ~ "Avian orthoavulavirus 9",
                                   Organism_Name == "Avian paramyxovirus 10" ~ "Avian metaavulavirus 10",
                                   Organism_Name == "Metaavulavirus peixense" ~ "Avian paramyxovirus 15",
                                   grepl("Mumps",Organism_Name,ignore.case = T) ~ "Mumps orthorubulavirus",
                                   grepl("Suis",Organism_Name,ignore.case = T) ~ "Orthorubulavirus suis",
                                   grepl("michoacan",Organism_Name,ignore.case = T) ~ "Orthorubulavirus suis",
                                   .default = Organism_Name))

write_xlsx(paramyxo_nuc_metadata,"paramyxovirus_metadata_cleaned_organism_names_newseqs.xlsx")

# read in old metadata, bind in new and write a new dataset containing both old and new sequence metadata
paramyxo_nuc_metadata_old <- read_xlsx("../../../3_host_assignment/paramyxovirus_metadata_cleaned_organism_names.xlsx")
paramyxo_nuc_metadata_full <- rbind(paramyxo_nuc_metadata,paramyxo_nuc_metadata_old)

write_xlsx(paramyxo_nuc_metadata_full,"paramyxovirus_metadata_cleaned_organism_names_full.xlsx")

# filter for only viral species not present in the original lookup table
new_metadata <- paramyxo_nuc_metadata %>% 
  filter(Organism_Name %nin% lookup_tbl$Organism_Name)

# function to get list of metadata hosts for every unique viral species
list_hosts <- function(x){
  target_virus <- new_metadata %>%
    filter(Organism_Name == x)

  return(unique(target_virus$Host))
}

virus_host_lists <- lapply(unique(new_metadata$Organism_Name),list_hosts)

no_host_vec <- vector()

# list the virus species for which no hosts are given in the metadata
for(i in 1:length(virus_host_lists)){
  if(F %in% is.na(virus_host_lists[[i]])){
  } else {
    no_host_vec <- c(no_host_vec,unique(new_metadata$Organism_Name)[i])
  }
}

# six new viral species present in the benchmark data
# all have clear host assignments in the metadata

# we now have a clean list of taxon names we can convert to IDs
nuc_host_taxa <- name2taxid(new_metadata$Host,db="ncbi",out_type="summary")

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

# this inelegant solution was implemented because the specific taxonomic ranks available varies widely across taxa
# no simple way to automate extraction particular taxonomic rank except explicitly by name
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
} else {print("No duplications")} 

# join this back to the metadata
colnames(nuc_host_taxa) <- c("Host","Host_id","Host_rank_class","Host_rank_order")

new_metadata <- left_join(new_metadata,nuc_host_taxa,by="Host")

# rm(nuc_missing,nuc_missing_host_taxa,nuc_missing_short,nuc_still_missing)

# add these new species to generate an expanded virus/host lookup table --------------------------------------
# all new species have clear host assignments
# note that the seal-associated virus may be dietary (in which case likely a fish-infecting virus)
bind_tbl <- new_metadata %>%
  select(Organism_Name,Host_rank_class,Host_rank_order) %>% 
  dplyr::rename("Host_Class" = Host_rank_class,"Host_Order" = Host_rank_order)

virus_host_lookup <- rbind(bind_tbl,lookup_tbl %>% select(-Host_label))

# assign a final host label: order level for mammalian viruses, class level for other viruses
virus_host_lookup %<>% 
  mutate(Host_label = case_when(Host_Class == "Aves" ~ "Aves",
                                Host_Class %in% c("Chondrichthyes","Actinopteri") ~ "Fish",
                                Host_Class == "Lepidosauria" ~ "Reptilia",
                                Host_Order == "Crocodylia" ~ "Reptilia", # this is done at order level as class assignment appears to have gone wrong - I confirmed that this virus was collected from caiman so this assignment is correct
                               .default = Host_Order
  ))

write_csv(virus_host_lookup,"virus_host_lookup_updated.csv")

# read and process nucleotide sequences ----------------------------------------------

# read in a TSV containing list of sequences for each cluster
clusters <- read.delim("../2_protein_orf_processing/0.9_0.4_clusters.tsv",header=F)

# for each cluster, get the accessions for every member
colnames(clusters) <- c("ref_accessions","members")

# assign organism names using the cleaned IDs in the metadata
virus_names <- paramyxo_nuc_metadata_full %>% 
  dplyr::select(Accession,Organism_Name) %>% 
  distinct()

# there are meant to be 2 joins here
clusters <- left_join(clusters,virus_names,join_by(ref_accessions == Accession))
clusters <- left_join(clusters,virus_names,join_by(members == Accession))

colnames(clusters) <- c("ref_accessions","members","ref_organism","member_organism")

host_names <- virus_host_lookup %>% 
  dplyr::select(Organism_Name,Host_label)

clusters <- left_join(clusters,host_names,join_by(ref_organism == Organism_Name))
clusters <- left_join(clusters,host_names,join_by(member_organism == Organism_Name))

colnames(clusters) <- c("ref_accessions","members","ref_organism","member_organism","ref_host_label","member_host_label")

# filter for any clusters containing viral species with different host labels
cluster_hosts <- function(x){
  
  target_cluster <- clusters %>%
    dplyr::filter(ref_accessions == x) %>% 
    drop_na(member_organism,ref_organism) # where an organism name is not available for specific sequences in the metadata, we ignore them
  
  return(unique(c(target_cluster$ref_host_label,target_cluster$member_host_label)))
}

host_cluster_list <- lapply(unique(clusters$ref_accessions),cluster_hosts)

multi_host_clusters <- vector()

for(i in 1:length(host_cluster_list)){
  if(F %in% is.na(host_cluster_list[[i]]) & length(host_cluster_list[[i]]) == T){
  } else {
    multi_host_clusters <- c(multi_host_clusters,unique(clusters$ref_accessions)[i])
  }
}

### no new multihost clusters - all new sequences are in non-conflicting clusters. Drop mammalian orthorubulavirus 5 as in previous dataset
# PV890885.1 and PZ154824.1 are both singleton clusters of mammalian orthorubulavirus 5

# make these corrections
clusters %<>% 
  mutate(Host_rank = ref_host_label)

# remove the generalist/indeterminate host clusters
clusters %<>% 
  filter(ref_accessions %nin% c("PV890885.1","PZ154824.1"))

clusters %<>% 
  distinct(ref_accessions,.keep_all = T) %>% 
  dplyr::select(ref_accessions,Host_rank)

write_csv(clusters,"all_cluster_hosts.csv")

# also remove these from the FASTA files
# whole genomes
whole_genome <- readDNAStringSet("../2_protein_orf_processing/centroid_0.9_0.4_seqs.fasta")

whole_genome <- whole_genome[!grepl(paste("PV890885.1","PZ154824.1",sep="|"),names(whole_genome))]

writeXStringSet(whole_genome,"centroid_0.9_0.4_seqs_clean.fasta")

# orf pseudogenome
pseudo_genome <- readDNAStringSet("../2_protein_orf_processing/complete_orf_pseudogenomes.fasta")

pseudo_genome <- pseudo_genome[!grepl(paste("PV890885.1","PZ154824.1",sep="|"),names(pseudo_genome))]

writeXStringSet(pseudo_genome,"complete_orf_pseudogenomes_clean.fasta")

# individual orfs
orf_seqs <- readDNAStringSet("../2_protein_orf_processing/nucleocapsid_orfs.fasta")

orf_seqs <- orf_seqs[!grepl(paste("PV890885.1","PZ154824.1",sep="|"),names(orf_seqs))]

writeXStringSet(orf_seqs,"nucleocapsid_orfs_clean.fasta")

orf_seqs <- readDNAStringSet("../2_protein_orf_processing/matrix_orfs.fasta")

orf_seqs <- orf_seqs[!grepl(paste("PV890885.1","PZ154824.1",sep="|"),names(orf_seqs))]

writeXStringSet(orf_seqs,"matrix_orfs_clean.fasta")

orf_seqs <- readDNAStringSet("../2_protein_orf_processing/fusion_orfs.fasta")

orf_seqs <- orf_seqs[!grepl(paste("PV890885.1","PZ154824.1",sep="|"),names(orf_seqs))]

writeXStringSet(orf_seqs,"fusion_orfs_clean.fasta")

orf_seqs <- readDNAStringSet("../2_protein_orf_processing/attachment_orfs.fasta")

orf_seqs <- orf_seqs[!grepl(paste("PV890885.1","PZ154824.1",sep="|"),names(orf_seqs))]

writeXStringSet(orf_seqs,"attachment_orfs_clean.fasta")

orf_seqs <- readDNAStringSet("../2_protein_orf_processing/polymerase_orfs.fasta")

orf_seqs <- orf_seqs[!grepl(paste("PV890885.1","PZ154824.1",sep="|"),names(orf_seqs))]

writeXStringSet(orf_seqs,"polymerase_orfs_clean.fasta")

# also remove from UTRs

# orf_seqs <- readDNAStringSet("../2_protein_orf_processing/3prime_leader_UTRs.fasta")
# 
# orf_seqs <- orf_seqs[!grepl(paste("PV890885.1","PZ154824.1",sep="|"),names(orf_seqs))]
# 
# writeXStringSet(orf_seqs,"3prime_leader_UTRs_clean.fasta")
# 
# orf_seqs <- readDNAStringSet("../2_protein_orf_processing/5prime_trailer_UTRs.fasta")
# 
# orf_seqs <- orf_seqs[!grepl(paste("PV890885.1","PZ154824.1",sep="|"),names(orf_seqs))]
# 
# writeXStringSet(orf_seqs,"5prime_trailer_UTRs_clean.fasta")
