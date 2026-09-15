library(tidyverse)
library(writexl)
library(readxl)
library(ORFik)
library(Biostrings)
library(GenomicFeatures)
library(magrittr)

`%nin%` <- Negate(`%in%`)

# IMPORTANT NOTE #
# all ORF assignments and validations were conducted for a specific clustering
# assignments were first conducted using strict (granular) clustering parameters to maximise the number of cluster centroids represented
# however, clustering with different parameters may still result in different cluster centroids that may require further manual assignment
# this is not avoidable - manual checks are essential to generate a high-quality dataset

# Extract ORFs ------------------------------------------------------------

cluster_centroids <- readDNAStringSet("./centroid_0.9_0.4_seqs.fasta")

# find ORFs - resolve fuzzy codons that can be unambiguously converted to stops
orfs <- findORFsFasta(cluster_centroids,startCodon = "ATG",stopCodon = paste(stopDefinition(transl_table = 1),"TRA","TAR",sep="|"))

# drop negative strand orfs
orfs <- orfs[strandBool(orfs)] 

long_orfs <- orfs[readWidths(orfs) > 400]

TEST_orfs <- orfs %>% as.data.frame()


# finding lots of ORFs
long_orfs_df <- long_orfs %>% as.data.frame()

# check for any sequences which are not returning ORFs
names(cluster_centroids)[str_split_i(names(cluster_centroids),"\\s\\|",i = 1) %nin% unique(long_orfs_df$seqnames) == T] 

long_orfs_df %>% group_by(seqnames) %>% tally %>% pull(n) %>% table
long_orfs_df %>% group_by(seqnames) %>% tally %>% arrange(-n)

# read in metadata table including the host information for all paramyxoviruses
paramyxo_nuc_metadata <- read_csv("../1_nucleotide_processing/1_initial_cleaning/1_quick_clean/11158_paramyxoviridae_complete_metadata_05112025_17062026.csv")

paramyxo_nuc_metadata <- paramyxo_nuc_metadata %>%
  dplyr::select(Accession,Organism_Name)

# join in virus names
long_orfs_df <- left_join(long_orfs_df,paramyxo_nuc_metadata,join_by("seqnames" == "Accession"))

# Processing validation refseqs -------------------------------------------

# we will download all available paramyxovirus refseqs and use these to validate our detected ORFs where they are available
prot_refseqs <- readAAStringSet("./paramyxoviridae_prot_refseqs_170626.fasta")

refseq_ids <- str_split_i(names(prot_refseqs),pattern = "\\s\\|",i=1)

refseq_ids[refseq_ids %in% unique(long_orfs_df$seqnames)] # check if we have any of these refseqs in our data

refseq_names <- str_split_i(names(prot_refseqs),pattern = "\\s\\|",i=2)

refseq_prot_names <- tibble(Protein = str_split_i(refseq_names,pattern = "\\s\\[",i=1))
refseq_species_names <- str_extract(refseq_names,pattern="\\[.*\\]")
refseq_species_names <- str_sub(refseq_species_names,start=2,end=-2)

refseq_names <- tibble(ref_organism = refseq_species_names)

refseq_names$ref_organism <- stringr::str_to_sentence(refseq_names$ref_organism)

# do some cleaning of organism names
refseq_names %<>%
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
                                  .default = ref_organism))

refseq_species_names <- refseq_names$ref_organism

# clean up protein names
refseq_prot_names %<>% mutate(Protein = case_when(  
  grepl("polymerase|RdRp|large|^L\\s|^L$|\\sL$", Protein, ignore.case = TRUE) ~ "Polymerase",
  grepl("matrix|^M\\s|^M$|\\sM$|\\sM\\s|matirx", Protein, ignore.case = TRUE) ~ "Matrix",
  grepl("nucleo|^N\\s|^N$|\\sN$|\\sN\\s|nucelo|nucle|nuleo|^NP\\s|^NP$", Protein, ignore.case = TRUE) ~ "Nucleocapsid",
  # we are combining all attachment protein classes (HN, H and G into one)
  grepl("^G\\s|^G$|\\sG$|hemagglutinin|haemagglutinin|hemagglutanin|hemagglutintin|hemaggluttinin|hemaggutinin|hemaglutinin|hemagluttinin|hemagglutin|neuraminidase|^HN\\s|^HN$|\\sHN\\s|\\sHN$|^H\\s|^H$|\\sH$|attachment", Protein, ignore.case = TRUE) ~ "Attachment",
  grepl("fusion|^F\\s|^F$|\\sF$|\\sF\\s", Protein, ignore.case = TRUE) ~ "Fusion",
  grepl("hydrophobic|^SH\\s|hydrofobic|hydropobic|^SH$", Protein, ignore.case = TRUE) ~ "Hydrophobic protein",
  # grepl("hemagglutinin|neuraminidase|^HN\\s|^HN$|\\sHN\\s", Protein, ignore.case = TRUE) ~ "HN", 
  grepl("phospho|^P\\s|^P$|\\sP$|^P\\'|^P\\-|\\sP\\s|phos", Protein, ignore.case = TRUE) ~ "Phosphoprotein",
  grepl("membrane|receptor|^TM\\s|^TM$|\\sTM\\s|\\sTM$", Protein, ignore.case = TRUE) ~ "Transmembrane protein",
  # grepl("^C\\s|^C$|\\sC$|^C\\'|^C\\-|\\sC\\s", Protein, ignore.case = TRUE) ~ "C protein",
  grepl("^V\\s|^V$|\\sV$|^V\\-|\\sV\\s|^W\\s|^W$|\\sW$|^W\\-|\\sW\\s|V\\/W|V\\/C|V\\/P|^C\\s|^C$|\\sC$|^C\\'|^C\\-|\\sC\\s|accessory|^D\\s|^I\\s|^Y1\\s|^Y2\\s|\\sY1\\s|\\sY2\\s", Protein, ignore.case = TRUE) ~ "Accessory",
)) %>% replace_na(list(Protein = "other"))

# paste back the full names to the sequences
names(prot_refseqs) <- paste(refseq_ids,refseq_prot_names$Protein,refseq_species_names,sep = "|")

# filter for only canonical proteins
# canonical_prots <- c("Polymerase","Matrix","Nucleocapsid","Attachment","Fusion","Phosphoprotein")
canonical_prots <- c("Polymerase","Matrix","Nucleocapsid","Attachment","Fusion")

prot_refseqs <- prot_refseqs[grepl(paste(canonical_prots,collapse="|"),names(prot_refseqs))]

# we now have a list of all ref canonical protein seqs for our selected validation viruses
# more than 1 set of refseqs are available for some viruses - remove duplicated sequences

# make a tibble for easy processing (we only need one refseq for each species)
refseqs_df <- tibble(data.frame(prot_refseqs))

refseqs_df$ID <- str_split_i(names(prot_refseqs),pattern = "\\|",i=1)
refseqs_df$Protein <- str_split_i(names(prot_refseqs),pattern = "\\|",i=2)
refseqs_df$Species <- str_split_i(names(prot_refseqs),pattern = "\\|",i=3)

# pull the source nucleotides from the metadata
refseqs_metadata <- read_csv("paramyxoviridae_prot_refseq_metadata_170626.csv")

refseqs_metadata %<>% 
  dplyr::select(Accession,Nucleotide)

refseqs_df <- left_join(refseqs_df,refseqs_metadata,by = join_by(ID == Accession))

# this pipe arranges by group size to ensure that wherever possible refseqs for each species come from the same nucleotide deposition
refseqs_df %<>% 
  group_by(Nucleotide) %>% 
  mutate(completeness = n()) %>% 
  arrange(desc(completeness))

# this pipe selects only unique protein/species combinations so we have a single refseq for each target protein for each viral species
refseqs_df <- refseqs_df %>% 
  group_by(Species,Protein) %>% 
  mutate(group_id = cur_group_id()) %>% 
  distinct(group_id,.keep_all = T)

# check if we have a full set of proteins for each species
refseqs_df %>% ungroup() %>% dplyr::count(Species) %>% arrange(n) # there should be 5 entries for each species
# we have all 5 proteins for most species for which any refseqs are available
# a few only have sequences for 4 proteins

# reassign prot_refseqs stringset using the unique sequences
prot_refseqs <- AAStringSet(x = refseqs_df$prot_refseqs)
names(prot_refseqs) <- paste(refseqs_df$ID,refseqs_df$Protein,refseqs_df$Species,sep = "|")

# Assigning ORFs ---------------------------------------------------------
centroid_names <- names(cluster_centroids)
# rename cluster centroids with only ID for easy matching
names(cluster_centroids) <- str_split_i(names(cluster_centroids),pattern = "\\s\\|",i=1)

# extract ORFs and give them a label based on length and location in genome

# approximate canonical relative lengths of proteins (from viral zone) for reference
nucleocapsid <- (525*3)/15100
matrix <- (335*3)/15100
fusion <- (550*3)/15100
attachment <- (617*3)/15100
polymerase <- (2183*3)/15100

long_orfs_df %<>% 
  group_by(seqnames) %>% 
  arrange(start,.by_group = T)

long_orfs_df %<>%
  group_by(seqnames) %>% 
  mutate(seq_length = width(cluster_centroids[seqnames])) %>% 
  mutate(rel_start=start/seq_length,rel_end=end/seq_length,rel_width=width/seq_length)

# NOTE #
# code block below was used during validation to avoid re-validating previously assessed sequences
# we used new_clusts instead of long_orfs_df during the validation
########################
# # read in list of sequences validated for the 0.8_0.8 clustering
# clusters_8_8 <- readAAStringSet("/Users/jamieherzig/Documents/supervised_training/paramyxovirus/complete_genomes/complete_orf_pseudogenomes.fasta")
# clusters_8_8 <- names(clusters_8_8)
# 
# 
# new_clusts <- long_orfs_df %>% 
#   filter(seqnames %nin% clusters_8_8)
# 
# # filter for sequences >17,000 as possible jeilongviruses (this will also include some henipas and related which have large genomes)
# view(new_clusts %>% 
#        filter(seq_length > 17000))
# 
# # this is probably all the new jeilongs (don't think any will be shorter than 17k bp but not certain)
# new_jeilong <- c("OQ379183.1","OL409126.1","KY370098.1","JN689227.1","AY900001.1")
########################

# split dataset into jeilong-type viruses and others
jeilong_set <- long_orfs_df %>% 
  ungroup() %>% 
  filter(grepl("jeilong|Jeilong|Beilong",Organism_Name) == T | 
           grepl("memana",Organism_Name) == T)

big_set <- long_orfs_df %>% 
  ungroup() %>% 
  filter(grepl(paste(jeilong_set$seqnames,collapse="|"), seqnames) == F)

# assign genes for jeilong-type viruses
jeilong_set %<>%
  dplyr::mutate(Protein = case_when(rel_start < 0.05 & rel_width > 0.07 & rel_width < 0.13 ~ "Nucleocapsid",
                                    rel_start > 0.15 & rel_end < 0.3 & rel_width > 0.04 & rel_width < 0.09 ~ "Matrix",
                                    rel_start > 0.22 & rel_end < 0.4 & rel_width > 0.07 & rel_width < 0.17 ~ "Fusion",
                                    rel_start > 0.37 & rel_end < 0.7 & rel_width > 0.09 & rel_width < 0.28 ~ "Attachment",
                                    rel_start > 0.55 & rel_width > 0.28 & rel_width < 0.5 ~ "Polymerase",
                                    .default="Undefined"))

# no manual corrections needed for this set of jeilongviruses

# assign genes for other paramyxoviruses
big_set %<>%
  ungroup() %>% 
  dplyr::mutate(Protein = case_when(rel_start < 0.1 & rel_width > 0.07 & rel_width < 0.13 ~ "Nucleocapsid",
                                    rel_start > 0.17 & rel_end < 0.35 & rel_width > 0.04 & rel_width < 0.09 ~ "Matrix",
                                    rel_start > 0.22 & rel_end < 0.5 & rel_width > 0.07 & rel_width < 0.17 ~ "Fusion",
                                    rel_start > 0.37 & rel_end < 0.62 & rel_width > 0.09 & rel_width < 0.17 ~ "Attachment",
                                    rel_start > 0.45 & rel_width > 0.32 & rel_width < 0.5 ~ "Polymerase",
                                    .default="Undefined"))

# manual corrections
big_set %<>% 
  dplyr::mutate(Protein = case_when(seqnames == "PX657369.1" & start == 4461 ~ "Undefined", # two matrix proteins assigned, correct assignment determined by pairwise alignment to ref
                                    seqnames == "PX682252.1" & start == 1239 ~ "Undefined", # two nucleocapsid assigned, correct assignment determined by pairwise alignment to bat paramyxovirus virus ref
                                    # PX841364.1 is a highly divergent seal-associated virus, challenging to validate
                                    seqnames == "PX841364.1" & start == 1449 ~ "Undefined", # two nucleocapsids assigned, no clear alignment to known seal or fish viruses but putative assignment made based on genome architecture (two ORFs of similar length in this region typically nucleocapsid followed by phosphoprotein)
                                    seqnames == "PX841364.1" & start == 4152 ~ "Undefined", # two matrix proteins assigned, no clear alignment to known seal or fish viruses but putative assignment made based on genome architecture (two ORFs of similar length in this region often seen with first, longer ORF validated as matrix in other viruses. Length of this ORF also more consistent with other paramyxovirus matrix proteins)
                                    .default = Protein))


# Note #
# codeblock below allows for alignment and inspection of individual ORFs against a reference
###############################
### WORKING ALIGNMENT
# current_seq <- cluster_centroids["PX841364.1"]
# 
# target_orfs <- big_set %>%
#   filter(seqnames == "PX841364.1")
# 
# orfs <- mapply(subseq,current_seq,start=target_orfs$start,end=target_orfs$end)
# 
# names(orfs) <- paste(target_orfs$Protein,"target",sep = " | ")
# 
# orfs <- DNAStringSet(orfs)
# 
# prot_orfs <- translate(orfs,if.fuzzy.codon = "solve")
# 
# alignment <- pwalign::pairwiseAlignment(prot_refseqs[grepl("Achimota pararubulavirus 3",names(prot_refseqs)) & grepl("Matrix",names(prot_refseqs))],prot_orfs[5])
# 
# alignment <- c(pwalign::alignedPattern(alignment),pwalign::alignedSubject(alignment))
# 
# DECIPHER::BrowseSeqs(alignment)
###############################

# bind back together
long_orfs_df <- rbind(jeilong_set,big_set)

long_orfs_df %<>%
  filter(Protein != "Undefined")

rm(orfs,long_orfs,jeilong_set,big_set)

# Validating ORFs ---------------------------------------------------------
# function to extract all ORFs from validation set genomes
extract_validation_orfs <- function(seq_id){
  
  current_seq <- cluster_centroids[seq_id]
  
  target_orfs <- long_orfs_df %>% 
    filter(seqnames == names(current_seq))
  
  orfs <- mapply(subseq,current_seq,start=target_orfs$start,end=target_orfs$end)
  
  names(orfs) <- paste(target_orfs$Protein,seq_id,sep = "|")
  
  orfs <- DNAStringSet(orfs)
  
  prot_orfs <- translate(orfs,if.fuzzy.codon = "solve")
  
  # we have a list of orfs, now align to the reference
  target_refseqs <- prot_refseqs[grepl(paste0(unique(target_orfs$Organism_Name),"$"),names(prot_refseqs))] # note the dollar sign anchor to avoid catching other numbered viral species
  
  # temporarily name the sequences by protein
  prot_orfs_prots <- str_split_i(names(prot_orfs),pattern = "\\|",i=1)
  prot_orfs_id <- str_split_i(names(prot_orfs),pattern = "\\|",i=2)

  target_id <- str_split_i(names(target_refseqs),pattern = "\\|",i=1)
  target_prots <- str_split_i(names(target_refseqs),pattern = "\\|",i=2)
  target_prots_species <- str_split_i(names(target_refseqs),pattern = "\\|",i=3)

  target_names_df <- tibble("id" = paste(target_prots_species,target_id,sep = "|"),"prots" = target_prots)
  
  names(prot_orfs) <- prot_orfs_prots
  names(target_refseqs) <- target_prots
  
  # align by extracting the proteins from target reference and subject with the same name
  # in some cases a reference is missing and will return an empty object, hence try wrappers
  try(N_alignment <- pwalign::pairwiseAlignment(target_refseqs[grepl(names(prot_orfs[1]),names(target_refseqs))],prot_orfs[1]),silent = T)
  try(M_alignment <- pwalign::pairwiseAlignment(target_refseqs[grepl(names(prot_orfs[2]),names(target_refseqs))],prot_orfs[2]),silent = T)
  try(F_alignment <- pwalign::pairwiseAlignment(target_refseqs[grepl(names(prot_orfs[3]),names(target_refseqs))],prot_orfs[3]),silent = T)
  try(HN_alignment <- pwalign::pairwiseAlignment(target_refseqs[grepl(names(prot_orfs[4]),names(target_refseqs))],prot_orfs[4]),silent = T)
  try(L_alignment <- pwalign::pairwiseAlignment(target_refseqs[grepl(names(prot_orfs[5]),names(target_refseqs))],prot_orfs[5]),silent = T)
  
  scores <- c(score(N_alignment),score(M_alignment),score(F_alignment),score(HN_alignment),score(L_alignment))
  
  # generate aligned stringsets for inspection with DECIPHER
  try(N_alignment <- c(pwalign::alignedPattern(N_alignment),pwalign::alignedSubject(N_alignment)),silent = T)
  try(M_alignment <- c(pwalign::alignedPattern(M_alignment),pwalign::alignedSubject(M_alignment)),silent = T)
  try(F_alignment <- c(pwalign::alignedPattern(F_alignment),pwalign::alignedSubject(F_alignment)),silent = T)
  try(HN_alignment <- c(pwalign::alignedPattern(HN_alignment),pwalign::alignedSubject(HN_alignment)),silent = T)
  try(L_alignment <- c(pwalign::alignedPattern(L_alignment),pwalign::alignedSubject(L_alignment)),silent = T)
  
  # rename with full IDs
  try(names(N_alignment) <- c(paste(target_names_df %>% filter(prots == names(N_alignment[1])) %>% dplyr::select(id)),paste(prot_orfs_prots[1],prot_orfs_id[1],sep="|")),silent = T)
  try(names(M_alignment) <- c(paste(target_names_df %>% filter(prots == names(M_alignment[1])) %>% dplyr::select(id)),paste(prot_orfs_prots[2],prot_orfs_id[2],sep="|")),silent = T)
  try(names(F_alignment) <- c(paste(target_names_df %>% filter(prots == names(F_alignment[1])) %>% dplyr::select(id)),paste(prot_orfs_prots[3],prot_orfs_id[3],sep="|")),silent = T)
  try(names(HN_alignment) <- c(paste(target_names_df %>% filter(prots == names(HN_alignment[1])) %>% dplyr::select(id)),paste(prot_orfs_prots[4],prot_orfs_id[4],sep="|")),silent = T)
  try(names(L_alignment) <- c(paste(target_names_df %>% filter(prots == names(L_alignment[1])) %>% dplyr::select(id)),paste(prot_orfs_prots[5],prot_orfs_id[5],sep="|")),silent = T)
  
  alignment_views <- list(N_alignment,M_alignment,F_alignment,HN_alignment,L_alignment)
  
  # return the orfs, a list of alignment scores, and the alignment view objects for manual validation if necessary
  return(list(orfs,scores,alignment_views))  
}
# 
### TESTING
# seq_id <- unique(validation_set$seqnames)[5]
# seq_id <- "FJ986192.2"

# generate a validation set containing only the viruses we are validating against refseqs
# repeat the cleaning of organism names we did for the refseqs
long_orfs_df$Organism_Name <- stringr::str_to_sentence(long_orfs_df$Organism_Name)

long_orfs_df %<>%
  mutate(Organism_Name = case_when(Organism_Name == "Paraavulavirus hongkongense" ~ "Avian paramyxovirus 4",
                                  Organism_Name == "Metaavulavirus hongkongense" ~ "Avian metaavulavirus 6",
                                  Organism_Name == "Avian paramyxovirus 6" ~ "Avian metaavulavirus 6",
                                  Organism_Name == "Avian paramyxovirus penguin/falkland islands/324/2007" ~ "Avian paramyxovirus 10",
                                  Organism_Name == "Jeilongvirus beilongi" ~ "Beilong virus",
                                  grepl("Rinderpest",Organism_Name,ignore.case = T) ~ "Rinderpest morbillivirus",
                                  Organism_Name == "Canine distemper virus" ~ "Morbillivirus canis",
                                  Organism_Name == "Human parainfluenza virus 4a" ~ "Human orthorubulavirus 4",
                                  Organism_Name == "Human parainfluenza virus 4b" ~ "Human orthorubulavirus 4",
                                  Organism_Name == "Small ruminant morbillivirus" ~ "Morbillivirus caprinae",
                                  grepl("Peste",Organism_Name,ignore.case = T) ~ "Morbillivirus caprinae",
                                  Organism_Name == "Avian paramyxovirus 1" ~ "Avian orthoavulavirus 1",
                                  Organism_Name == "Orthoavulavirus javaense" ~ "Avian orthoavulavirus 1",
                                  Organism_Name == "Pigeon paramyxovirus 1" ~ "Avian orthoavulavirus 1",
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
                                  Organism_Name == "Orthoavulavirus newyorkense" ~ "Avian orthoavulavirus 9",
                                  Organism_Name == "Avian paramyxovirus 10" ~ "Avian metaavulavirus 10",
                                  Organism_Name == "Metaavulavirus peixense" ~ "Avian paramyxovirus 15",
                                  grepl("Mumps",Organism_Name,ignore.case = T) ~ "Mumps orthorubulavirus",
                                  .default = Organism_Name))

validation_set <- long_orfs_df %>% 
  filter(Organism_Name %in% refseq_species_names)

print(length(unique(validation_set$Organism_Name))) # 71 unique viruses for which we have a refseq available
print(length(unique(validation_set$seqnames))) # 134 unique sequences for which we have a refseq available

### NOTE ###
# this function is coded to handle missing ORFs in the reference sequences
# however, it is hard coded to assume every subject being validated has 5 ORFs and will throw errors if this is not the case

validation_alignments <- lapply(unique(validation_set$seqnames),extract_validation_orfs)

# summarise results
# scores
# set up an empty tibble to bind to
scoring_tbl <- tibble("Scores"=NA,"Protein"=NA,"Reference_ID"=NA,"Subject_ID"=NA,"Organism"=NA,.rows=0)

for(i in 1:length(validation_alignments)){
  
  subject_ids <- c(try(str_split_i(names(validation_alignments[[i]][[3]][[1]][1]),pattern = "\\|",i=2),silent = T),
               try(str_split_i(names(validation_alignments[[i]][[3]][[2]][1]),pattern = "\\|",i=2),silent = T),
               try(str_split_i(names(validation_alignments[[i]][[3]][[3]][1]),pattern = "\\|",i=2),silent = T),
               try(str_split_i(names(validation_alignments[[i]][[3]][[4]][1]),pattern = "\\|",i=2),silent = T),
               try(str_split_i(names(validation_alignments[[i]][[3]][[5]][1]),pattern = "\\|",i=2),silent = T))
  
  subject_ids <- subject_ids[!grepl("Error",subject_ids)]
  
  organism_names <- c(try(str_split_i(names(validation_alignments[[i]][[3]][[1]][1]),pattern = "\\|",i=1),silent = T),
                      try(str_split_i(names(validation_alignments[[i]][[3]][[2]][1]),pattern = "\\|",i=1),silent = T),
                      try(str_split_i(names(validation_alignments[[i]][[3]][[3]][1]),pattern = "\\|",i=1),silent = T),
                      try(str_split_i(names(validation_alignments[[i]][[3]][[4]][1]),pattern = "\\|",i=1),silent = T),
                      try(str_split_i(names(validation_alignments[[i]][[3]][[5]][1]),pattern = "\\|",i=1),silent = T))
  
  organism_names <- organism_names[!grepl("Error",organism_names)]
  
  protein_names <- c(try(str_split_i(names(validation_alignments[[i]][[3]][[1]][2]),pattern = "\\|",i=1),silent = T),
                     try(str_split_i(names(validation_alignments[[i]][[3]][[2]][2]),pattern = "\\|",i=1),silent = T),
                     try(str_split_i(names(validation_alignments[[i]][[3]][[3]][2]),pattern = "\\|",i=1),silent = T),
                     try(str_split_i(names(validation_alignments[[i]][[3]][[4]][2]),pattern = "\\|",i=1),silent = T),
                     try(str_split_i(names(validation_alignments[[i]][[3]][[5]][2]),pattern = "\\|",i=1),silent = T))
  
  protein_names <- protein_names[!grepl("Error",protein_names)]
  
  ref_ids <- c(try(str_split_i(names(validation_alignments[[i]][[3]][[1]][2]),pattern = "\\|",i=2),silent = T),
                     try(str_split_i(names(validation_alignments[[i]][[3]][[2]][2]),pattern = "\\|",i=2),silent = T),
                     try(str_split_i(names(validation_alignments[[i]][[3]][[3]][2]),pattern = "\\|",i=2),silent = T),
                     try(str_split_i(names(validation_alignments[[i]][[3]][[4]][2]),pattern = "\\|",i=2),silent = T),
                     try(str_split_i(names(validation_alignments[[i]][[3]][[5]][2]),pattern = "\\|",i=2),silent = T))
  
  ref_ids <- ref_ids[!grepl("Error",ref_ids)]
  
  temp_tbl <- tibble("Scores"=validation_alignments[[i]][[2]],"Protein"=protein_names,"Reference_ID"=ref_ids,"Subject_ID"=subject_ids,"Organism"=organism_names,"Index"=i)
  
  scoring_tbl <- rbind(scoring_tbl,temp_tbl)
}

# inspect any with negative scores
neg_scores <- scoring_tbl %>% 
  filter(Scores <= 0)

# no negative scores for this set, all ref alignments are good

# Write outputs -----------------------------------------------------------
# function to extract all ORFs from each protein
extract_orfs <- function(seq_id){
  
  current_seq <- cluster_centroids[seq_id]
  
  target_orfs <- long_orfs_df %>% 
    filter(seqnames == names(current_seq))
  
  orfs <- mapply(subseq,current_seq,start=target_orfs$start,end=target_orfs$end)
  
  orfs <- DNAStringSet(orfs)
  
  names(orfs) <- paste(seq_id,target_orfs$Protein,sep="|")
  
  # return individual ORFs (for mafft alignment)
  return(orfs)
}

# function to extract all ORFs from each protein
extract_all_orfs <- function(seq_id){
  
  current_seq <- cluster_centroids[seq_id]
  
  target_orfs <- long_orfs_df %>% 
    filter(seqnames == names(current_seq))
  
  orfs <- mapply(subseq,current_seq,start=target_orfs$start,end=target_orfs$end)
  
  orfs <- DNAStringSet(orfs)
  
  # make a new DNAstring object containing the 5 ORFs concatenated together as a pseudogenome
  # try calls handle cases where not all canonical ORFs are present
  all_orfs <- paste0(try(as.character(orfs[1])),
                     try(as.character(orfs[2])),
                     try(as.character(orfs[3])),
                     try(as.character(orfs[4])),
                     try(as.character(orfs[5])))
  
  # remove any non-canonical characters
  all_orfs <- paste(str_extract_all(all_orfs,c("A|T|G|C"))[[1]],collapse="")
  
  all_orfs <- DNAStringSet(all_orfs)
  
  names(all_orfs) <- seq_id
  
  # return the pseudogenome
  return(all_orfs)
}

orf_indiv <- do.call(c,lapply(unique(long_orfs_df$seqnames),extract_orfs))
orf_genomes <- do.call(c,lapply(unique(long_orfs_df$seqnames),extract_all_orfs))

# split into individual protein stringsets
orf_N <- orf_indiv[grepl("Nucleocapsid",names(orf_indiv))]
orf_M <- orf_indiv[grepl("Matrix",names(orf_indiv))]
orf_F <- orf_indiv[grepl("Fusion",names(orf_indiv))]
orf_HN <- orf_indiv[grepl("Attachment",names(orf_indiv))]
orf_L <- orf_indiv[grepl("Polymerase",names(orf_indiv))]

writeXStringSet(orf_N,"nucleocapsid_orfs.fasta")
writeXStringSet(orf_M,"matrix_orfs.fasta")
writeXStringSet(orf_F,"fusion_orfs.fasta")
writeXStringSet(orf_HN,"attachment_orfs.fasta")
writeXStringSet(orf_L,"polymerase_orfs.fasta")

writeXStringSet(orf_genomes,"complete_orf_pseudogenomes.fasta")

# write the long_orfs_df table so we can locate these orfs within the genome in future
write_xlsx(long_orfs_df,"all_orfs_table.xlsx")
