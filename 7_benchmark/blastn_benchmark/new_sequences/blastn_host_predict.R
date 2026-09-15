library(tidyverse)
library(readxl)
library(magrittr)
library(e1071)
library(igraph)
library(rcartocolor)

### this script extracts the top hits following blast of ORFs
`%nin%` <- Negate(`%in%`)

N_orf_hits <- read_tsv("Norfs_dcmegablast_results.txt",col_names = F)
M_orf_hits <- read_tsv("Morfs_dcmegablast_results.txt",col_names = F)
F_orf_hits <- read_tsv("Forfs_dcmegablast_results.txt",col_names = F)
HN_orf_hits <- read_tsv("HNorfs_dcmegablast_results.txt",col_names = F)
L_orf_hits <- read_tsv("Lorfs_dcmegablast_results.txt",col_names = F)

col_names <- c("query_id","target_id","shared_seqid","alignment_length","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bitscore")

colnames(N_orf_hits) <- col_names
colnames(M_orf_hits) <- col_names
colnames(F_orf_hits) <- col_names
colnames(HN_orf_hits) <- col_names
colnames(L_orf_hits) <- col_names

# filter for best hits for each query based on bit score and evalue
find_top_hits <- function(x){
  
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

# read in virus-host lookup table and clean organism name metadata
lookup_tbl <- read_csv("../../../7_benchmark/new_sequences/3_host_assignment/virus_host_lookup_updated.csv")
organism_names <- read_xlsx("../../../7_benchmark/new_sequences/3_host_assignment/paramyxovirus_metadata_cleaned_organism_names_full.xlsx")

lookup_tbl %<>%
  select(Organism_Name,Host_label)

organism_names %<>%
  select(Accession,Organism_Name)

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

colnames(final_preds) <- c("Accession","Predicted_label")

results_pred <- inner_join(results_pred,final_preds)

results_pred %<>%
  mutate(Outcome = case_when(Host_label == Predicted_label ~ "Correct",
                            .default = "Incorrect"))

write_csv(results_pred,"new_seqs_blast_pred_results.csv")
# calculate accuracy
accuracy <- results_pred %>% 
  count(Outcome) 

accuracy <- accuracy$n[2]/sum(accuracy$n)

table(results_pred$Outcome)

# 
# # plot outputs on graph ---------------------------------------------------
# # read in graph data
# edgelist <- read_delim("/Users/jamieherzig/Documents/clean_run/5_phylogeny/3_generate_graph/graph_data/patristic_paraymxo_graph.edgelist",col_names = F)
# edgelist <- edgelist+1
# 
# g <- graph_from_edgelist(as.matrix(edgelist),directed = F)
# 
# # read in metadata
# branch_weights <- read_csv("/Users/jamieherzig/Documents/clean_run/5_phylogeny/3_generate_graph/graph_data/patristic_paramyxo_graph_branch_weights.csv")
# colnames(branch_weights) <- c("Node","Weight")
# 
# node_names <- read_csv("/Users/jamieherzig/Documents/clean_run/5_phylogeny/3_generate_graph/graph_data/patristic_paramyxo_graph_node_names.csv")
# colnames(node_names) <- c("Node","Accession")
# 
# E(g)$weight <- branch_weights$Weight
# V(g)$names <- node_names$Accession
# 
# # hosts
# host_labels <- read_csv("/Users/jamieherzig/Documents/clean_run/3_host_assignment/all_cluster_hosts.csv")
# 
# # read in layout
# layout_tbl <- read_csv("/Users/jamieherzig/Documents/clean_run/5_phylogeny/3_generate_graph/graph_phylo_plots/final_manual_layout.csv")
# 
# layout <- layout_tbl %>% select(x_coord,y_coord)
# layout <- as.matrix(layout)
# 
# # use a join to order the prediction results correctly
# node_res <- left_join(node_names,results_pred)
# 
# node_res %<>%
#   mutate(result = case_when(result == T ~ "Correct",
#                             .default = "Incorrect"))
# 
# V(g)$outcome <- node_res$result
# 
# node_res <- left_join(node_res,host_labels,join_by("Accession" == "ref_accessions"))
# 
# node_res$Host_rank <- factor(node_res$Host_rank,levels=c("Primates","Rodentia","Carnivora","Chiroptera","Artiodactyla","Perissodactyla","Scandentia","Pholidota","Eulipotyphla",
#                                                                                "Aves","Reptilia","Fish","Undetermined"))
# V(g)$label <- node_res$Host_rank
# 
# # set up a colour palette for plotting
# # palette <- c(carto_pal(12, "Safe"),"grey1")
# palette <- carto_pal(12, "Safe")
# 
# mycol = palette[node_res$Host_rank]
# 
# shapes <- c("circle","square")
# 
# myshape <- shapes[factor(node_res$result,levels = c("Correct","Incorrect"))]
# 
# png(filename="blast_results_graph.png",width=2000,heigh=4000,units="px")
# plot(g, vertex.size=3, vertex.label=NA,vertex.color=mycol,layout = layout,edge.width=3)
# legend('topright',legend=levels(factor((V(g)$label))),col = palette,pch = 16,cex=4)
# legend('bottomright',legend = c("Correct","Incorrect"),pch = c(1,0),cex=4)
# dev.off()
# 
# png(filename="blast_results_graph_labelled.png",width=2000,heigh=4000,units="px")
# plot(g, vertex.size=3, vertex.label=V(g)$names,vertex.color=mycol,vertex.shape=myshape,layout = layout,edge.width=3)
# legend('topright',legend=levels(factor((V(g)$label))),col = palette,pch = 16,cex=4)
# legend('bottomright',legend = c("Correct","Incorrect"),pch = c(1,0),cex=4)
# dev.off()
# 
# view(node_res %>% filter(result == "Incorrect"))
