library(tidyverse)
library(readxl)
library(viridis)
library(magrittr)
library(igraph)
library(rcartocolor)

# define a function to apply function to every combination of 2 arguments
cmapply <- function(FUN, ..., MoreArgs = NULL){  
  # expand a grid of all argument combinations  
  l <- expand.grid(..., stringsAsFactors=FALSE)    
  # apply the function  
  .mapply(FUN=FUN, dots=unname(l), MoreArgs = MoreArgs)}

best_tune <- read_csv("../../UTR_kfold_cv/simple_feats/results_summary/best_tune_acc_allUTR.csv")

read_best_tunes <- function(feat,condition,nhidden,nlayers){
  blocked_res <- read_csv(paste0("./feature_sets/results_",condition,"_",feat,"/loo_results_",nhidden,"_",nlayers,"convlayers_weighted.csv"))
  blocked_res$feat <- feat
  blocked_res$condition <- condition
  
  return(blocked_res)
}

best_tune_results <- do.call(rbind,mapply(read_best_tunes,
                                          feat = best_tune$feat,
                                          condition = best_tune$condition,
                                          nhidden = best_tune$nhidden,
                                          nlayers = best_tune$nlayers,
                                          SIMPLIFY = F))

best_tune_results %<>% 
  mutate(Correct = case_when(Observed == Predicted ~ T,
                             .default = F))

view(best_tune_results %>% 
  count(feat,condition,Correct))

# read in UTR lengths to see if they correlate with results
UTR_lengths <- read_csv("../../../../4_feature_generation/3_genomic_feats/all_UTR_lengths.csv")

leader_res <- best_tune_results %>%
  filter(condition == "leader")

leader_res <- left_join(leader_res,UTR_lengths)

cor(x = leader_res$leader_length,y = leader_res$Correct,method = "pearson")

trailer_res <- best_tune_results %>%
  filter(condition == "trailer")

trailer_res <- left_join(trailer_res,UTR_lengths)

cor(x = trailer_res$trailer_length,y = trailer_res$Correct,method = "spearman")

view(leader_res %>% 
       filter(Correct == F) %>% 
       select(Observed,Predicted,feat,leader_length))

view(trailer_res %>% 
       filter(Correct == F) %>% 
       select(Observed,Predicted,feat,trailer_length))

view(best_tune_results %>% 
  filter(Correct == F) %>% 
  select(Observed,Predicted,feat))

test <- leader_res |> filter(condition == "leader",feat == "kmer6only",Correct == F)

test <- test$leader_length

wilcox.test(UTR_lengths$leader_length,test)

