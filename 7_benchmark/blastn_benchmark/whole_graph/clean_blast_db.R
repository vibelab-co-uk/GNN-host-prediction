library(tidyverse)
library(readxl)
library(magrittr)

# quick script to identify all sequences for which we do not have a host
organism_names <- read_xlsx("/Users/jamieherzig/Documents/clean_run/3_host_assignment/paramyxovirus_metadata_cleaned_organism_names.xlsx")

# we remove 
# Mammalian orthorubulavirus 5
# Phyllostomus bat morbillivirus
# Marmoset morbillivirus

drop_seqs <- organism_names %>%
  filter(Organism_Name %in% c("Mammalian orthorubulavirus 5","Phyllostomus bat morbillivirus","Marmoset morbillivirus")) %>% 
  select(Accession)

write_lines(drop_seqs$Accession,"seqs_to_drop.txt")

