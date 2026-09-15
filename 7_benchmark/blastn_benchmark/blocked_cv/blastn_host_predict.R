library(tidyverse)
library(readxl)
library(magrittr)
library(e1071)
library(igraph)
library(rcartocolor)

`%nin%` <- Negate(`%in%`)

find_top_hits <- function(x){
  
  if(nrow(x) == 0){
    return()
  }
  
  x %<>% 
    separate(query_id,into=c("query_id","orf"),sep="\\|") %>%  # split to give a separate orf id and query id
    filter(query_id != target_id) # drop self matches
  
  # filter for top 1% bit score & lowest 1% evalue
  top_hits <- x %>% 
    group_by(query_id) %>% 
    filter(bitscore >= quantile(bitscore,probs = 0.99) & evalue <= quantile(evalue,probs = 0.01)) # using a quantile cut off results in variable numbers of matches
  # this is realistic as the number and closeness of matches will be variable 
  
  return(top_hits)
}

### this script extracts the top hits following blast of ORFs
analyse_blast <- function(x){
  
  N_orf_hits <- read_tsv(paste0(x,"/Norfs_dcmegablast_results.txt"),col_names = F)
  M_orf_hits <- read_tsv(paste0(x,"/Morfs_dcmegablast_results.txt"),col_names = F)
  F_orf_hits <- read_tsv(paste0(x,"/Forfs_dcmegablast_results.txt"),col_names = F)
  HN_orf_hits <- read_tsv(paste0(x,"/HNorfs_dcmegablast_results.txt"),col_names = F)
  L_orf_hits <- read_tsv(paste0(x,"/Lorfs_dcmegablast_results.txt"),col_names = F)
  
  col_names <- c("query_id","target_id","shared_seqid","alignment_length","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bitscore")
  
  colnames(N_orf_hits) <- col_names
  colnames(M_orf_hits) <- col_names
  colnames(F_orf_hits) <- col_names
  colnames(HN_orf_hits) <- col_names
  colnames(L_orf_hits) <- col_names
  
  # filter for best hits for each query based on bit score and evalue
  all_hits <- do.call(rbind,lapply(list(N_orf_hits,M_orf_hits,F_orf_hits,HN_orf_hits,L_orf_hits),find_top_hits))
  
  # calculate the skew in the bitscore
  all_hits %<>% 
    group_by(query_id,orf) %>% 
    mutate(skew = skewness(bitscore))
  
  # for highly skewed data, perform another round of quantile filtering to ensure we only use closest matches
  skew_quant <- quantile(all_hits$skew,na.rm=T,probs = seq(from=0.1,to = 1,by = 0.1))
  
  low_skew <- all_hits %>% 
    filter(skew < 0.5 | is.na(skew) == T) # skew is NA where only 1 hit is found, keep these records
  
  # for highly skewed data, we perform another round of quantile filtering depending on the severity of the skew
  high_skew <- all_hits %>% 
    group_by(query_id,orf) %>% 
    filter(case_when(skew >= 0.5 & skew < 1 ~ bitscore >= quantile(bitscore,probs = 0.75)[[1]],
                     skew >= 1 & skew < 2 ~ bitscore >= quantile(bitscore,probs = 0.9)[[1]],
                     skew >= 2 ~ bitscore >= quantile(bitscore,probs = 0.95)[[1]])) 
  
  all_hits <- rbind(low_skew,high_skew)
  
  # bind into the blast hits
  all_hits <- inner_join(all_hits,organism_names, join_by(target_id == Accession))
  
  all_hits <- inner_join(all_hits,lookup_tbl)
  
  # remove viruses for which we could not determine a host label
  all_hits %<>%
    filter(is.na(Host_label) != T)
  
  # calculate modal host from near hits across all 5 ORFs -------------------
  # set up a results table
  results_pred <- tibble("Accession" = unique(all_hits$query_id))
  results_pred <- inner_join(results_pred,organism_names)
  results_pred <- inner_join(results_pred,lookup_tbl)
  
  orf_preds <- all_hits %>% 
    count(query_id,orf,Host_label) %>% 
    filter(n == max(n)) %>% 
    ungroup()
  
  final_preds <- orf_preds %>% 
    count(query_id,Host_label) %>% 
    group_by(query_id) %>% 
    filter(n == max(n))
  
  # for cases where we have a tied vote, we predict the host with highest sum bitscore
  tied_votes <- final_preds %>% 
    count(query_id) %>% 
    filter(n > 1)
  
  if(nrow(tied_votes)>0){
    tied_votes <- all_hits %>% 
      filter(query_id %in% tied_votes$query_id)
    
    tied_votes %<>% 
      group_by(query_id,Host_label) %>% 
      summarise(bitscore_sum = sum(bitscore)) %>% 
      filter(bitscore_sum == max(bitscore_sum)) %>% 
      select(-bitscore_sum)
    
    # remove extra columns, drop tied queries and rebind in resolved prediction
    final_preds %<>%
      select(-n) %>% 
      filter(query_id %nin% tied_votes$query_id)
    
    final_preds <- rbind(final_preds,tied_votes)
  } else {
    final_preds %<>%
      select(-n)
  }
  
  colnames(final_preds) <- c("Accession","Predicted_label")
  
  results_pred <- inner_join(results_pred,final_preds)
  
  results_pred %<>%
    mutate(result = case_when(Host_label == Predicted_label ~ T,
                              .default = F))
  

  # write results
  write_csv(results_pred,paste0(x,"/results.csv"))
  results_pred$Subgraph <- str_extract(x,"\\d+")
  return(results_pred)
}

# read in virus-host lookup table and clean organism name metadata
lookup_tbl <- read_csv("../../../3_host_assignment/virus_host_lookup.csv")
organism_names <- read_xlsx("../../../3_host_assignment/paramyxovirus_metadata_cleaned_organism_names.xlsx")

lookup_tbl %<>%
  select(Organism_Name,Host_label)

organism_names %<>%
  select(Accession,Organism_Name)

# apply function
directories <- list.files(path = ".",pattern = "subgraph*")

all_res <- do.call(rbind,lapply(directories,analyse_blast))

subgraph_acc <- all_res %>% 
  count(Subgraph,result)

subgraph_acc %<>% 
  group_by(Subgraph) %>% 
  mutate(total = sum(n)) %>% 
  pivot_wider(names_from = result, values_from = n)

colnames(subgraph_acc) <- c( "Subgraph","total","Incorrect","Correct" )

# calculate accuracy for all runs
subgraph_acc %<>% 
  mutate(Correct = replace_na(Correct, 0),Incorrect = replace_na(Incorrect, 0)) %>% 
  mutate(Accuracy = Correct/total) %>% 
  select(Subgraph,Accuracy)

write_csv(subgraph_acc,"blast_blocked_accuracy_metrics.csv")
