library(Biostrings)
library(magrittr)

paramyxo_genomes <- readDNAStringSet("11158_paramyxoviridae_complete_nucseqs_21112025.fasta")

# filter out sequences with >1% ambiguous positions
paramyxo_genomes <- paramyxo_genomes[letterFrequency(paramyxo_genomes, letters = "N") < ceiling(width(paramyxo_genomes)*0.01)]

# check quantiles to find sensible threshold to remove fragmentary sequences
quants <- quantile(width(paramyxo_genomes),probs = seq(0,1,0.05))

# 20% looks sensible
paramyxo_genomes <- paramyxo_genomes[width(paramyxo_genomes) >= quants[5]]

# remove frameshifted measles genomes and transgenic strains
drop_seqs <- paramyxo_genomes[grepl("Measles",names(paramyxo_genomes)) & 
                           grepl("UNVERIFIED",names(paramyxo_genomes)) &
                           grepl("Marseille",names(paramyxo_genomes)) | # catches all low quality/fragmented sequences f
                             grepl("GFP",names(paramyxo_genomes)) |
                             grepl("LC187306.1",names(paramyxo_genomes)) | # transgenic vaccine strains
                             grepl("MH144178.1",names(paramyxo_genomes)) | # transgenic vaccine strains
                             grepl("MW090971.1",names(paramyxo_genomes)) | # transgenic vaccine strains
                             grepl("OP456084.1",names(paramyxo_genomes))] # low quality/fragmented sequence

# write dropped seqs for reference
drop_seqs_ref <- tibble("seqs" = names(drop_seqs))

drop_seqs_ref %<>%
  separate_wider_delim(seqs,delim = "|",names = c("Accession","Sequence description"))

write_csv(drop_seqs_ref,"dropped_sequences.csv")

`%nin%` <- Negate(`%in%`)

paramyxo_genomes <- paramyxo_genomes[names(paramyxo_genomes) %nin% names(drop_seqs)]

# we'll also remove sequences deposited after 05/11/2025 to ensure consistency with previous analyses
old_genomes <- readDNAStringSet("./old_seqs/11158_paramyxoviridae_complete_nucseqs_05112025.fasta")
paramyxo_genomes <- paramyxo_genomes[names(paramyxo_genomes) %in% names(old_genomes)]

# write clean set
writeXStringSet(paramyxo_genomes,"11158_paramyxoviridae_cleaned_nucseqs_21112025.fasta")
