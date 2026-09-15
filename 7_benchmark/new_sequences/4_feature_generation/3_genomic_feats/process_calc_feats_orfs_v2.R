library(tidyverse)
library(Biostrings)
library(coRdon)
library(protr)
library(magrittr)
library(writexl)
library(readxl)
library(fishualize)

source("./feat_calc_funcs.R")
`%nin%` <- Negate(`%in%`)

# this script writes individual feature sets that can then be combined to perform macro feature selection

# read clustered paramyxovirus sequences and calculate genomic features
paramyxo_nuc <- readDNAStringSet("../../3_host_assignment/centroid_0.9_0.4_seqs_clean.fasta")

id_list <- str_split_i(names(paramyxo_nuc),pattern = "\\s\\|",i=1)

names(paramyxo_nuc) <- paste(id_list,"Genomic",sep="|")

# read in lists of unaligned ORFs
N_orfs <- readDNAStringSet("../../3_host_assignment/nucleocapsid_orfs_clean.fasta")
M_orfs <- readDNAStringSet("../../3_host_assignment/matrix_orfs_clean.fasta")
F_orfs <- readDNAStringSet("../../3_host_assignment/fusion_orfs_clean.fasta")
HN_orfs <- readDNAStringSet("../../3_host_assignment/attachment_orfs_clean.fasta")
L_orfs <- readDNAStringSet("../../3_host_assignment/polymerase_orfs_clean.fasta")

all_orfs <- DNAStringSet(c(paramyxo_nuc,N_orfs,M_orfs,F_orfs,HN_orfs,L_orfs))

# calculate genomic features
setwd("../../4_feature_generation/3_genomic_feats/nucleotide_feats")

# calculate nucleotide feats for all orfs as well as sequence level
# calc and write gc content
output <- calc_composition_prop_gc(all_orfs)

output %<>%
  separate_wider_delim(nuc_id,names=c("Accession","Protein"),delim = "|") %>% 
  pivot_wider(names_from=Protein,values_from=-c(Accession,Protein))
  
write_csv(output,"paramyxo_gc_prop_features.csv")

# calc and write random pseudofeature sets with the same dimensions as gc content feature
sample_set <- seq(from=30,to=70,by=0.001) # generate a random set to sample

for(i in 1:5){
  rand_set <- output # assign the gc content feature table
  
  # replace values with random numbers, maintaining the same dimension size
  rand_set$Genomic <-  sample(sample_set,size=nrow(output),replace=TRUE) 
  rand_set$Nucleocapsid <-  sample(sample_set,size=nrow(output),replace=TRUE) 
  rand_set$Matrix <-  sample(sample_set,size=nrow(output),replace=TRUE) 
  rand_set$Fusion <-  sample(sample_set,size=nrow(output),replace=TRUE) 
  rand_set$Attachment <-  sample(sample_set,size=nrow(output),replace=TRUE) 
  rand_set$Polymerase <-  sample(sample_set,size=nrow(output),replace=TRUE) 
  
  # write the random feature set
  write_csv(rand_set,paste0("paramyxo_rand_features_set",i,".csv"))
}

# calc and write nucleotide content
output <- calc_composition_prop_nt(all_orfs)

output %<>%
  separate_wider_delim(nuc_id,names=c("Accession","Protein"),delim = "|") %>% 
  pivot_wider(names_from=Protein,values_from=-c(Accession,Protein))

write_csv(output,"paramyxo_nt_prop_features.csv")

# calc and write dinucleotide content
output <- calc_composition_prop_dint(all_orfs)

output %<>%
  separate_wider_delim(nuc_id,names=c("Accession","Protein"),delim = "|") %>% 
  pivot_wider(names_from=Protein,values_from=-c(Accession,Protein))

write_csv(output,"paramyxo_dint_prop_features.csv")

# calc and write 3mer content
output <- calc_composition_prop_3mer(all_orfs)

output %<>%
  separate_wider_delim(nuc_id,names=c("Accession","Protein"),delim = "|") %>% 
  pivot_wider(names_from=Protein,values_from=-c(Accession,Protein))

write_csv(output,"paramyxo_3mer_prop_features.csv")

# calc and write 4mer content
output <- calc_composition_prop_4mer(all_orfs)

output %<>%
  separate_wider_delim(nuc_id,names=c("Accession","Protein"),delim = "|") %>% 
  pivot_wider(names_from=Protein,values_from=-c(Accession,Protein))

write_csv(output,"paramyxo_4mer_prop_features.csv")

# calc and write 5mer content
output <- calc_composition_prop_5mer(all_orfs)

output %<>%
  separate_wider_delim(nuc_id,names=c("Accession","Protein"),delim = "|") %>% 
  pivot_wider(names_from=Protein,values_from=-c(Accession,Protein))

write_csv(output,"paramyxo_5mer_prop_features.csv")

# calc and write 6mer content
output <- calc_composition_prop_6mer(all_orfs)

output %<>%
  separate_wider_delim(nuc_id,names=c("Accession","Protein"),delim = "|") %>% 
  pivot_wider(names_from=Protein,values_from=-c(Accession,Protein))

write_csv(output,"paramyxo_6mer_prop_features.csv")

setwd("..")

############################################################################# 
############################################################################# 

# read and process amino acid metadata & sequences -----------------------------------------------------
# read in protein orfs from codon-corrected data in protr format
all_orfs <- readFASTA("../1_ambiguous_codon_assignment/all_orfs_ambig_corrected.fasta")

# remove the stop codons (protr needs sequences with only 20 generic AAs)
all_orfs_names <- names(all_orfs)

all_orfs <- as.list(str_replace_all(all_orfs,pattern="\\*",""))

names(all_orfs) <- all_orfs_names

# ========================================================
aa_features <- calc_composition_prop_aac(all_orfs)

aa_features %<>%
  separate_wider_delim(prot_id,names=c("Accession","Protein"),delim = "|") %>% 
  pivot_wider(names_from=Protein,values_from=-c(Accession,Protein)) %>% 
  replace(is.na(.), 0)

write_csv(aa_features,"./protein_feats/paramyxo_aac_prop_features.csv")

# for some of the count features, we replace NAs with 0s where a particular count does not occur

aa_features <- calc_composition_prop_dp(all_orfs)

aa_features %<>%
  separate_wider_delim(prot_id,names=c("Accession","Protein"),delim = "|") %>% 
  pivot_wider(names_from=Protein,values_from=-c(Accession,Protein)) %>% 
  replace(is.na(.), 0)

write_csv(aa_features,"./protein_feats/paramyxo_dp_prop_features.csv")


aa_features <- calc_moranac_lag30(all_orfs)

aa_features %<>%
  separate_wider_delim(prot_id,names=c("Accession","Protein"),delim = "|") %>% 
  pivot_wider(names_from=Protein,values_from=-c(Accession,Protein)) %>% 
  replace(is.na(.), 0)

write_csv(aa_features,"./protein_feats/paramyxo_moranac_lag30_features.csv")

aa_features <- calc_expanded_moranac_lag30(all_orfs)

aa_features %<>%
  separate_wider_delim(prot_id,names=c("Accession","Protein"),delim = "|") %>% 
  pivot_wider(names_from=Protein,values_from=-c(Accession,Protein)) %>% 
  replace(is.na(.), 0)

write_csv(aa_features,"./protein_feats/paramyxo_expanded_moranac_lag30_features.csv")


aa_features <- calc_moranac_lag100(all_orfs)

aa_features %<>%
  separate_wider_delim(prot_id,names=c("Accession","Protein"),delim = "|") %>% 
  pivot_wider(names_from=Protein,values_from=-c(Accession,Protein)) %>% 
  replace(is.na(.), 0)

write_csv(aa_features,"./protein_feats/paramyxo_moranac_lag100_features.csv")


aa_features <- calc_ctdc(all_orfs)

aa_features %<>%
  separate_wider_delim(prot_id,names=c("Accession","Protein"),delim = "|") %>% 
  pivot_wider(names_from=Protein,values_from=-c(Accession,Protein)) %>% 
  replace(is.na(.), 0)

write_csv(aa_features,"./protein_feats/paramyxo_ctdc_features.csv")

aa_features <- calc_ctdt(all_orfs)

aa_features %<>%
  separate_wider_delim(prot_id,names=c("Accession","Protein"),delim = "|") %>% 
  pivot_wider(names_from=Protein,values_from=-c(Accession,Protein)) %>% 
  replace(is.na(.), 0)

write_csv(aa_features,"./protein_feats/paramyxo_ctdt_features.csv")

aa_features <- calc_ctdd(all_orfs)

aa_features %<>%
  separate_wider_delim(prot_id,names=c("Accession","Protein"),delim = "|") %>% 
  pivot_wider(names_from=Protein,values_from=-c(Accession,Protein)) %>% 
  replace(is.na(.), 0)

write_csv(aa_features,"./protein_feats/paramyxo_ctdd_features.csv")

aa_features <- calc_ctriad(all_orfs)

aa_features %<>%
  separate_wider_delim(prot_id,names=c("Accession","Protein"),delim = "|") %>% 
  pivot_wider(names_from=Protein,values_from=-c(Accession,Protein)) %>% 
  replace(is.na(.), 0)

write_csv(aa_features,"./protein_feats/paramyxo_ctriad_features.csv")
