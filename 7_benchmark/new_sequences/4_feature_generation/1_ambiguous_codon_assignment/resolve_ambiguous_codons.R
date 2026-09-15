library(tidyverse)
library(Biostrings)
library(taxizedb)
library(protr)
library(magrittr)
library(writexl)
library(readxl)

# script to resolve ambiguous codons in cluster centroid sequences where possible
source("./ambiguous_codon_funcs.R")
`%nin%` <- Negate(`%in%`)

# read in lists of unaligned ORFs
N_orfs <- readDNAStringSet("../../3_host_assignment/nucleocapsid_orfs_clean.fasta")
M_orfs <- readDNAStringSet("../../3_host_assignment/matrix_orfs_clean.fasta")
F_orfs <- readDNAStringSet("../../3_host_assignment/fusion_orfs_clean.fasta")
HN_orfs <- readDNAStringSet("../../3_host_assignment/attachment_orfs_clean.fasta")
L_orfs <- readDNAStringSet("../../3_host_assignment/polymerase_orfs_clean.fasta")

# translate 
N_orfs <- translate(N_orfs,if.fuzzy.codon = "solve")
M_orfs <- translate(M_orfs,if.fuzzy.codon = "solve")
F_orfs <- translate(F_orfs,if.fuzzy.codon = "solve")
HN_orfs <- translate(HN_orfs,if.fuzzy.codon = "solve")
L_orfs <- translate(L_orfs,if.fuzzy.codon = "solve")

# rename with proteins in names
# bind all together
# afterwards add accession and into datatable for binding

all_orfs <- AAStringSet(c(N_orfs,M_orfs,F_orfs,HN_orfs,L_orfs))

writeXStringSet(all_orfs,"translated_orfs.fasta")

# read back in in protr format
all_orfs <- readFASTA("translated_orfs.fasta")

# remove the stop codons (protr needs sequences with only 20 generic AAs)
all_orfs_names <- names(all_orfs)

all_orfs <- as.list(str_replace_all(all_orfs,pattern="\\*",""))

# locate ambiguous sites
ambig_sites <- str_locate_all(all_orfs,pattern="X")

collect_ambig_sites <- function(x){
  tibble(start = ambig_sites[[x]][,1], end = ambig_sites[[x]][,2], ref_id = all_orfs_names[x])
}

# collect all ambiguous sites so we can later remove them from the ESM datatable
ambig_sites <- do.call(rbind,lapply(1:length(ambig_sites),collect_ambig_sites))

# for sequences with identified ambiguous sites, read in any other sequences from the cluster
# align these sequences, generate a consensus matrix and then replace the ambiguous AA with the consensus AA
unique(ambig_sites$ref_id)

orf_table <- read_xlsx("../../2_protein_orf_processing/all_orfs_table.xlsx")
ambig_sites$nuc_start <- ambig_sites$start * 3
ambig_sites$nuc_end <- ambig_sites$end * 3 + 2

### HERE HERE HERE 
# only reading in the new clusters duuuuh
# ... but that's fine we only need to do the new ones

# read in clustering table
clusters <- read.delim("../../2_protein_orf_processing/0.9_0.4_clusters.tsv",header=F)
colnames(clusters) <- c("ref_accessions","members")

# read in the full set of sequences
ambig_fasta <- readDNAStringSet("../../1_nucleotide_processing/1_initial_cleaning/2_combined_dataset_cleaning/11158_paramyxoviridae_cleaned_nucseqs_combined_nodups.fasta")

# filter for only clusters where centroids have ambiguous codons
ambig_seqs <- unique(str_split_i(ambig_sites$ref_id,"\\|",i = 1))

clusters %<>%
  filter(ref_accessions %in% ambig_seqs)

# filter the fasta for only sequences in these clusters
ambig_fasta <- ambig_fasta[grepl(paste(clusters$members,collapse = "|"),names(ambig_fasta))]

corrected_ambig_tbl <- do.call(rbind,lapply(ambig_seqs,resolve_ambig))

if(is.null(corrected_ambig_tbl) == F){

# if any of the determined consensus codons code for a stop, add them to a separate tbl for removal and drop them from the corrected tbl
drop_me <- corrected_ambig_tbl %>%
  filter(consensus_aa == "*")

corrected_ambig_tbl %<>%
  filter(consensus_aa != "*")

# for orfs where we couldn't find a consensus, drop ambiguous positions
names(all_orfs) <- all_orfs_names

# drop ambiguous amino acids from sequences present in ambig_seqs but not present in the corrected_ambig_tbl
all_orfs[grepl(paste(ambig_seqs[ambig_seqs %nin% corrected_ambig_tbl$seqnames],collapse="|"),names(all_orfs))] <- as.list(str_replace_all(all_orfs[grepl(paste(ambig_seqs[ambig_seqs %nin% corrected_ambig_tbl$seqnames],collapse="|"),names(all_orfs))],pattern="X",""))

# do the same for any where the consensus was a stop codon
for(i in 1:nrow(drop_me)){
  all_orfs[grepl(paste0(drop_me$seqnames[i],"|",drop_me$Protein[i]),names(all_orfs),fixed=T)] <- as.list(str_replace_all(all_orfs[grepl(paste0(drop_me$seqnames[i],"|",drop_me$Protein[i]),names(all_orfs),fixed=T)],pattern="X",""))
}

# replace the correct amino acids in the sequences we could find consensus sequences for
for(i in 1:nrow(corrected_ambig_tbl)){
  str_sub(all_orfs[grepl(corrected_ambig_tbl$ref_id[i],names(all_orfs),fixed=T)],corrected_ambig_tbl$start[i],corrected_ambig_tbl$end[i]) <- corrected_ambig_tbl$consensus_aa[i]
}

} else { # if we couldn't find any consensus, simply drop the ambiguous amino acids

  # drop ambiguous amino acids from sequences present in ambig_seqs but not present in the corrected_ambig_tbl
  all_orfs[grepl(paste(ambig_seqs[ambig_seqs %nin% corrected_ambig_tbl$seqnames],collapse="|"),names(all_orfs))] <- as.list(str_replace_all(all_orfs[grepl(paste(ambig_seqs[ambig_seqs %nin% corrected_ambig_tbl$seqnames],collapse="|"),names(all_orfs))],pattern="X",""))
  
}
# finally convert to biostrings format again and write updated fasta files for downstream use
all_orfs <- AAStringSet(c(unlist(unname(all_orfs))))
names(all_orfs) <- all_orfs_names

# these corrected amino acid sequences will be used to calculate protein genomic and PLM features
writeXStringSet(all_orfs,"./all_orfs_ambig_corrected.fasta")
