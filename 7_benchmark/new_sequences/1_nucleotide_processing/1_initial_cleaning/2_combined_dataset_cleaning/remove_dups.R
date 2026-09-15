library(Biostrings)
library(dplyr)
### this script removes duplicated sequences from the combined old + new fasta file
# note that dups will be handled through clustering in the main pipeline, but here they are removed as new refseqs may exist in the old dataset as non-refseqs

dup_IDs <- read_tsv("duplicated_seqs.txt",col_names = F)

colnames(dup_IDs) <- c("N","sequences")

new_names <- paste0("sequence_",1:max(dup_IDs$N))

# one set to drop from the combined (it doesn't matter which we keep)
dup_IDs_combined <- dup_IDs %>% 
  separate_wider_delim(sequences,delim = ", ",names = new_names,cols_remove=T,too_few = "align_start") %>% 
  dplyr::select(-c(N,sequence_1)) %>%  # keep the sequences in the first column, drop all others
  pivot_longer(cols = everything(),names_to = NULL,values_to = "sequences") %>% 
  drop_na()

combined_seqs <- readDNAStringSet("11158_paramyxoviridae_cleaned_nucseqs_combined.fasta")

drop_seqs <- combined_seqs[grepl(paste(dup_IDs_combined$sequences,collapse = "|"),names(combined_seqs))]

`%nin%` <- Negate(`%in%`)

combined_seqs <- combined_seqs[names(combined_seqs) %nin% names(drop_seqs)]

writeXStringSet(combined_seqs,"11158_paramyxoviridae_cleaned_nucseqs_combined_nodups.fasta")

# for a new set we need to make sure we are dropping the sequences from the new set
# e.g. the duplicated sequences we drop need to be those actually present in the new, not their duplicates in the old
# e.g. we first need to filter the duplicated sequences to find those that are present in the new seqs
# then we can drop those

new_seqs <- readDNAStringSet("../1_quick_clean/11158_paramyxoviridae_cleaned_nucseqs_05112025_17062026.fasta")

new_seq_ids <- str_split_i(names(new_seqs),pattern = "\\s\\|",i =1)

dup_IDs_new <- dup_IDs %>% 
  separate_wider_delim(sequences,delim = ", ",names = new_names,cols_remove=T,too_few = "align_start")

dup_IDs_new$index <- 1:nrow(dup_IDs_new)

dup_IDs_new %<>%
  pivot_longer(cols = starts_with("sequence"),names_to = NULL,values_to = "sequences") %>% 
  drop_na() %>% 
  mutate(new_seq = case_when(sequences %in% (new_seq_ids) == T ~ T,
                             .default = F))

# filter for only mixed groups (e.g where sequences are duplicated between old and new sets)
dup_IDs_new %<>% 
  group_by(index) %>% 
  filter(length(unique(new_seq)) > 1)

# they are all single duplications - simply drop the duplicates in the new set
dup_IDs_new %<>% 
  filter(new_seq == T)

drop_seqs <- combined_seqs[grepl(paste(dup_IDs_new$sequences,collapse = "|"),names(combined_seqs))]

`%nin%` <- Negate(`%in%`)

new_seqs <- new_seqs[names(new_seqs) %nin% names(drop_seqs)]

writeXStringSet(new_seqs,"../1_quick_clean/11158_paramyxoviridae_cleaned_nucseqs_05112025_17062026_nodups.fasta")







