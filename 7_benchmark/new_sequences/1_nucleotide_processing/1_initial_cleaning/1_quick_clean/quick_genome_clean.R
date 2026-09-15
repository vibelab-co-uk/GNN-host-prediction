library(Biostrings)

paramyxo_genomes <- readDNAStringSet("11158_paramyxoviridae_complete_nucseqs_05112025_17062026.fasta")

# filter out sequences with >1% ambiguous positions
paramyxo_genomes <- paramyxo_genomes[letterFrequency(paramyxo_genomes, letters = "N") < ceiling(width(paramyxo_genomes)*0.01)]

# check quantiles to find sensible threshold to remove fragmentary sequences
quants <- quantile(width(paramyxo_genomes),probs = seq(0,1,0.05))

# remove a number of fragmentary measles genomes, all others look good
paramyxo_genomes <- paramyxo_genomes[width(paramyxo_genomes) > quants[2]]

drop_seqs <- paramyxo_genomes[grepl("PZ011079.1",names(paramyxo_genomes))]
`%nin%` <- Negate(`%in%`)

paramyxo_genomes <- paramyxo_genomes[names(paramyxo_genomes) %nin% names(drop_seqs)]

# write clean set
writeXStringSet(paramyxo_genomes,"11158_paramyxoviridae_cleaned_nucseqs_05112025_17062026.fasta")
