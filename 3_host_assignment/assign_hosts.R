library(tidyverse)
library(rentrez)
library(Biostrings)
library(taxizedb)
library(coRdon)
library(protr)
library(magrittr)
library(writexl)
library(fishualize)

# read and process nucleotide metadata -----------------------------------------------------

# read in metadata table including the host information for all paramyxoviruses
paramyxo_nuc_metadata <- read_csv("../1_nucleotide_processing/1_initial_cleaning/11158_paramyxoviridae_complete_metadata_21112025.csv")

# clean virus species names
paramyxo_nuc_metadata$Organism_Name <- stringr::str_to_sentence(paramyxo_nuc_metadata$Organism_Name)

paramyxo_nuc_metadata %<>%
  mutate(Organism_Name = case_when(Organism_Name == "Paraavulavirus hongkongense" ~ "Avian paramyxovirus 4",
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

write_xlsx(paramyxo_nuc_metadata,"paramyxovirus_metadata_cleaned_organism_names.xlsx")

# function to get list of metadata hosts for every unique viral species
list_hosts <- function(x){
  target_virus <- paramyxo_nuc_metadata %>%
    filter(Organism_Name == x)

  return(unique(target_virus$Host))
}

virus_host_lists <- lapply(unique(paramyxo_nuc_metadata$Organism_Name),list_hosts)

no_host_vec <- vector()

# list the virus species for which no hosts are given in the metadata
for(i in 1:length(virus_host_lists)){
  if(F %in% is.na(virus_host_lists[[i]])){
  } else {
    no_host_vec <- c(no_host_vec,unique(paramyxo_nuc_metadata$Organism_Name)[i])
  }
}

# manually assign some hosts that are missing from the metadata but are clear in the literature
paramyxo_nuc_metadata <- paramyxo_nuc_metadata %>% 
  mutate(Host = case_when(Organism_Name == "Orthorubulavirus simiae" ~ "Macaca fascicularis", # lab strain of SV41 originally isolated from cynomologus monkey kidneys
                          Organism_Name == "Paraavulavirus neophemae" ~ "Aves", # lab strain of PMV-3/parakeet/Netherlands/449/75 grown in chick embryo
                          Organism_Name == "Tupaia paramyxovirus" ~ "Tupaia belangeri", # tupaia paramyxovirus isolated from tree shrews
                          Organism_Name == "Fer-de-lance virus" ~ "Bothrops atrox", # Fer-de-Lance virus, host given as fer-de-lance although it infects other reptiles as I believe this is where lab strain originated
                          .default = Host)
  )

# we now have at least one host in the metadata for every viral species, we can drop NA hosts without losing any info
paramyxo_nuc_metadata_nohost <- paramyxo_nuc_metadata %>%
  drop_na(Host)

### NOTE
# code block below used to manually identify and correct any names that are not correctly parsed by name2taxid
# with the cleaning carried out above, we get all of them for the paramyxo dataset so not currently needed
###

#########################################################
#########################################################
# now we want to find the parent taxonomy at given level for every host
# nuc_host_taxa <- name2taxid(paramyxo_nuc_metadata_nohost$Host,db="ncbi",out_type="summary")
# 
# check for any missing
# length(unique(paramyxo_nuc_metadata$Host))-nrow(nuc_host_taxa)
# 
# # look at them
# nuc_missing <- paramyxo_nuc_metadata %>%
#   filter(!(Host %in% nuc_host_taxa$name))
# 
# # make a duplicate to modify
# nuc_missing_short <- nuc_missing
# 
# # first just drop species name, this should fix most of them
# nuc_missing_short$Host <- str_extract(nuc_missing$Host,"^.*(?=\\s)")
# 
# nuc_missing_host_taxa <- name2taxid(nuc_missing_short$Host,db="ncbi",out_type="summary")
# 
# # check if any still missing
# length(unique(nuc_missing_short$Host))-nrow(nuc_missing_host_taxa)
# 
# # inspect and correct manually
# nuc_still_missing <- nuc_missing_short %>%
#   filter(!(Host %in% nuc_missing_host_taxa$name))
# 
# # correct all these in the main datatable
# # first drop species rank
# paramyxo_nuc_metadata <- paramyxo_nuc_metadata %>%
#   mutate(Host = case_when((Host %in% nuc_missing$Host) == T ~ str_extract(Host,"^.*(?=\\s)"),
#                           .default = Host
#   ))
# 
# # then manual fixes
# paramyxo_nuc_metadata <- paramyxo_nuc_metadata %>% 
#   mutate(Host = case_when(Host == "Scoterepens" ~ "Scotorepens", # typo in metadata
#                           Host == "Lissonycteris" ~ "Myonycteris", # angolan fruit bat, NCBI taxonomy considers it a member of myonycteris
#                           .default = Host)
#   )
#########################################################
#########################################################

# we now have a clean list of taxon names we can convert to IDs
nuc_host_taxa <- name2taxid(paramyxo_nuc_metadata_nohost$Host,db="ncbi",out_type="summary")
rm(paramyxo_nuc_metadata_nohost)

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
}

# manually correct these
nuc_host_taxa <- nuc_host_taxa %>% 
  filter(rank_Class != "Magnoliopsida") # this is a plant tax rank which shares a name with a pangolin tax rank

# join this back to the metadata
colnames(nuc_host_taxa) <- c("Host","Host_id","Host_rank_class","Host_rank_order")

paramyxo_nuc_metadata <- left_join(paramyxo_nuc_metadata,nuc_host_taxa,by="Host")

# rm(nuc_missing,nuc_missing_host_taxa,nuc_missing_short,nuc_still_missing)

# generate a virus/host lookup table --------------------------------------
virus_host_lookup <- tibble(Organism_Name=unique(paramyxo_nuc_metadata$Organism_Name))

# add the host taxonomic class for every unique virus in the dataset
class_assign <- function(x){
  
  target_virus <- paramyxo_nuc_metadata %>% 
    filter(Organism_Name == x) 
  
  target_virus <- tibble(Organism_Name = NA, Host_Class = unique(target_virus$Host_rank_class))
  
  target_virus$Organism_Name <- x
  
  target_virus %<>%
    drop_na()
  
  return(target_virus)
}

class_hosts <- do.call(rbind,lapply(virus_host_lookup$Organism_Name,class_assign))

# filter for hosts with unanimous assignment
single_class_hosts <- class_hosts %>%
  group_by(Organism_Name) %>% 
  filter(n() == 1)

virus_host_lookup <- left_join(virus_host_lookup,single_class_hosts)

# filter for only viruses with multiple recorded hosts at the class level
multi_class_hosts <- class_hosts %>%
  group_by(Organism_Name) %>% 
  filter(n() > 1)

# determine a true reservoir host class for these species where possible
virus_host_lookup %<>% 
  mutate(Host_Class = case_when(Organism_Name == "Avian orthoavulavirus 1" ~ "Aves", # NDV spills over to mammals semi-frequently, evolutionary origin host class clearly aves
                                Organism_Name == "Mammalian orthorubulavirus 5" ~ "Mammalia", # single sequence isolated from a lizard, evolutionary origin host class clearly mammalia
                               .default = Host_Class
  ))

# repeat this process with host orders
# add the host taxonomic order for every unique virus in the dataset
order_assign <- function(x){
  
  target_virus <- paramyxo_nuc_metadata %>% 
    filter(Organism_Name == x) 
  
  target_virus <- tibble(Organism_Name = NA, Host_Order = unique(target_virus$Host_rank_order))
  
  target_virus$Organism_Name <- x
  
  target_virus %<>%
    drop_na()
  
  return(target_virus)
}

order_hosts <- do.call(rbind,lapply(virus_host_lookup$Organism_Name,order_assign))

# filter for hosts with unanimous assignment
single_order_hosts <- order_hosts %>%
  group_by(Organism_Name) %>% 
  filter(n() == 1)

virus_host_lookup <- left_join(virus_host_lookup,single_order_hosts)

# filter for only viruses with multiple recorded hosts at the order level
multi_order_hosts <- order_hosts %>%
  group_by(Organism_Name) %>% 
  filter(n() > 1)

# determine a true reservoir host order for these species where possible
virus_host_lookup %<>% 
  mutate(Host_Order = case_when(Organism_Name == "Sosuga virus" ~ "Chiroptera", # reservoir species is the egyptian fruit bat with spillovers into human population
                                Organism_Name == "Morbillivirus caprinae" ~ "Artiodactyla", # PdPRV sometimes spills over into dogs, reservoir host artiodactyla
                                Organism_Name == "Morbillivirus canis" ~ "Carnivora", # CTD sometimes spills over into humans or livestock, reservoir host is carnivora
                                Organism_Name == "Rinderpest morbillivirus" ~ "Artiodactyla", # rinderpest sometimes spills over into lagomorphs, reservoir host artiodactyla
                                Organism_Name == "Henipavirus nipahense" ~ "Chiroptera", # nipah virus spills over into several other mammalian orders, reservoir hosts are all bats
                                Organism_Name == "Henipavirus hendraense" ~ "Chiroptera", # similar to nipah, hendra spills over readily but reservoir is 4 bat species
                                Organism_Name == "Langya virus" ~ "Eulipotyphla", # as above, spillovers to multiple species but good evidence shrews are the reservoir host
                                Organism_Name == "Respirovirus muris" ~ "Rodentia", # Sendai virus only infects mice as far as I know, one sequence isolated from marmoset lung but no publication available
                                .default = Host_Order
  ))


# assign a final host label: order level for mammalian viruses, class level for other viruses
virus_host_lookup %<>% 
  mutate(Host_label = case_when(Host_Class == "Aves" ~ "Aves",
                                Host_Class %in% c("Chondrichthyes","Actinopteri") ~ "Fish",
                                Host_Class == "Lepidosauria" ~ "Reptilia",
                                Host_Order == "Crocodylia" ~ "Reptilia", # this is done at order level as class assignment appears to have gone wrong - I confirmed that this virus was collected from caiman so this assignment is correct
                               .default = Host_Order
  ))

write_csv(virus_host_lookup,"virus_host_lookup.csv")

# read and process nucleotide sequences ----------------------------------------------

# read in a TSV containing list of sequences for each cluster
clusters <- read.delim("/Users/jamieherzig/Documents/clean_run/2_protein_orf_processing/0.9_0.4_clusters.tsv",header=F)

# for each cluster, get the accessions for every member
colnames(clusters) <- c("ref_accessions","members")

# assign organism names using the cleaned IDs in the metadata
virus_names <- paramyxo_nuc_metadata %>% 
  dplyr::select(Accession,Organism_Name)

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
    filter(ref_accessions == x) %>% 
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

# four clusters with viruses infecting multiple host labels/NAs
# KX060176.1 is a mammalian orthorubulavirus 5 cluster (also contains one goose orthorubulavirus). This is a generalist, drop
# MK593539.1 is a singleton cluster of mammalian orthorubulavirus 5
# MN306039.1 is a human respirovirus cluster. It contains a single sequence derived from pangolin, which was found to be human respirovirus 3. See DOI: 10.1080/22221751.2022.2086071
# MZ312422.1 is a cluster of morbilliviruses collected from Brazil. The centroid is from bat but all other sequences are collected over several years from marmosets. Origin host not clear, drop this cluster

# make these corrections
# assign "MN306039.1" as human-infecting
clusters %<>% 
  mutate(Host_rank = case_when(ref_accessions == "MN306039.1" ~ "Primates",
                                .default = ref_host_label
  ))

`%nin%` <- Negate(`%in%`)

# remove the generalist/indeterminate host clusters
clusters %<>% 
  filter(ref_accessions %nin% c("KX060176.1","MK593539.1","MZ312422.1"))

clusters %<>% 
  distinct(ref_accessions,.keep_all = T) %>% 
  dplyr::select(ref_accessions,Host_rank)

write_csv(clusters,"all_cluster_hosts.csv")

# also remove these from the FASTA files
# whole genomes
whole_genome <- readDNAStringSet("../2_protein_orf_processing/centroid_0.9_0.4_seqs.fasta")

whole_genome <- whole_genome[!grepl(paste("KX060176.1","MK593539.1","MZ312422.1",sep="|"),names(whole_genome))]

writeXStringSet(whole_genome,"centroid_0.9_0.4_seqs_clean.fasta")

# orf pseudogenome
pseudo_genome <- readDNAStringSet("../2_protein_orf_processing/complete_orf_pseudogenomes.fasta")

pseudo_genome <- pseudo_genome[!grepl(paste("KX060176.1","MK593539.1","MZ312422.1",sep="|"),names(pseudo_genome))]

writeXStringSet(pseudo_genome,"complete_orf_pseudogenomes_clean.fasta")

# individual orfs
orf_seqs <- readDNAStringSet("../2_protein_orf_processing/nucleocapsid_orfs.fasta")

orf_seqs <- orf_seqs[!grepl(paste("KX060176.1","MK593539.1","MZ312422.1",sep="|"),names(orf_seqs))]

writeXStringSet(orf_seqs,"nucleocapsid_orfs_clean.fasta")

orf_seqs <- readDNAStringSet("../2_protein_orf_processing/matrix_orfs.fasta")

orf_seqs <- orf_seqs[!grepl(paste("KX060176.1","MK593539.1","MZ312422.1",sep="|"),names(orf_seqs))]

writeXStringSet(orf_seqs,"matrix_orfs_clean.fasta")

orf_seqs <- readDNAStringSet("../2_protein_orf_processing/fusion_orfs.fasta")

orf_seqs <- orf_seqs[!grepl(paste("KX060176.1","MK593539.1","MZ312422.1",sep="|"),names(orf_seqs))]

writeXStringSet(orf_seqs,"fusion_orfs_clean.fasta")

orf_seqs <- readDNAStringSet("../2_protein_orf_processing/attachment_orfs.fasta")

orf_seqs <- orf_seqs[!grepl(paste("KX060176.1","MK593539.1","MZ312422.1",sep="|"),names(orf_seqs))]

writeXStringSet(orf_seqs,"attachment_orfs_clean.fasta")

orf_seqs <- readDNAStringSet("../2_protein_orf_processing/polymerase_orfs.fasta")

orf_seqs <- orf_seqs[!grepl(paste("KX060176.1","MK593539.1","MZ312422.1",sep="|"),names(orf_seqs))]

writeXStringSet(orf_seqs,"polymerase_orfs_clean.fasta")

# also remove from UTRs

orf_seqs <- readDNAStringSet("../2_protein_orf_processing/3prime_leader_UTRs.fasta")

orf_seqs <- orf_seqs[!grepl(paste("KX060176.1","MK593539.1","MZ312422.1",sep="|"),names(orf_seqs))]

writeXStringSet(orf_seqs,"3prime_leader_UTRs_clean.fasta")

orf_seqs <- readDNAStringSet("../2_protein_orf_processing/5prime_trailer_UTRs.fasta")

orf_seqs <- orf_seqs[!grepl(paste("KX060176.1","MK593539.1","MZ312422.1",sep="|"),names(orf_seqs))]

writeXStringSet(orf_seqs,"5prime_trailer_UTRs_clean.fasta")
