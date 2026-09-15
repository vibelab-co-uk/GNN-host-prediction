library(tidyverse)
library(readxl)
library(viridis)
library(magrittr)
library(igraph)
library(rcartocolor)

### this is a utility script to pull the best tunes to generate plots of incorrect nodes for 5-fold CV

`%nin%` <- Negate(`%in%`)

# read best tunes for all features ----------------------------------------
best_tune <- read_csv("./results_summary/best_tune_acc.csv")

# edit best tune tibble to include all rand sets
rand_best_tune <- tibble(nhidden = best_tune$nhidden[best_tune$feat == "rand"],nlayers = best_tune$nlayers[best_tune$feat == "rand"],feat = c("rand1","rand2","rand3","rand4","rand5"),mean_acc = best_tune$mean_acc[best_tune$feat == "rand"])

best_tune %<>%
  filter(feat != "rand")

best_tune <- rbind(best_tune,rand_best_tune)

# define a function to apply function to every combination of 2 arguments
cmapply <- function(FUN, ..., MoreArgs = NULL){  
  # expand a grid of all argument combinations  
  l <- expand.grid(..., stringsAsFactors=FALSE)    
  # apply the function  
  .mapply(FUN=FUN, dots=unname(l), MoreArgs = MoreArgs)}

read_best_tunes <- function(feat,nhidden,nlayers){
  blocked_res <- read_csv(paste0("./feature_sets/results_",feat,"/abl_output_",nhidden,"_",nlayers,".csv"))
  blocked_res$feat <- feat
  return(blocked_res)
}

best_tune_results <- tibble()

for(i in 1:nrow(best_tune)){
  current_feat <- best_tune$feat[i]
  tmp_results <- cmapply(read_best_tunes,
                         nhidden = list(as.vector(best_tune %>% filter(feat == current_feat) %>% select(nhidden) %>% unname()))[[1]][[1]],
                         nlayers = list(as.vector(best_tune %>% filter(feat == current_feat) %>% select(nlayers) %>% unname()))[[1]][[1]],
                         MoreArgs = list(feat = current_feat))
  
  best_tune_results <- rbind(best_tune_results,tmp_results[[1]])
}

# process results ---------------------------------------------------------

### block below calculates an ensemble prediction from the 5 CV folds
# calculate final ensemble prediction
# # make final prediction across ensemble folds for each feature
# cv_preds <- best_tune_results %>%
#   group_by(Accession,feat) %>% 
#   summarise(across(starts_with("tax"),~mean(.x)))
# 
# # process random pseudofeatures
# rand_results <- cv_preds %>%
#   filter(feat %in% c("rand1","rand2","rand3","rand4","rand5"))
# 
# # take mean results across the five random pseudofeatures
# rand_results %<>%
#   group_by(Accession) %>%
#   summarise(across(starts_with("tax"),~mean(.x)))
# 
# rand_results$feat <- "rand"
# 
# feat_results <- cv_preds %>%
#   filter(feat %nin% c("rand1","rand2","rand3","rand4","rand5"))
# 
# cv_preds <- rbind(feat_results,rand_results)
# 
# # make final predictions
# cv_preds %<>% 
#   group_by(Accession,feat) %>% 
#   mutate(Prediction_score = max(across(starts_with("tax")))) %>% 
#   # subsequent mutate assigns the name of the tax column with the max score (thereby giving the prediction)
#   mutate(Predicted = str_sub(colnames(cv_preds)[grepl("^tax",colnames(cv_preds))][which.max(across(starts_with("tax")))],start = 5))
# 
# cv_preds$feat <- factor(cv_preds$feat,levels = c("gconly","dintonly","kmer3only","aaconly","dponly","rand"))
# 
# # bind back in some metadata
# cv_preds <- left_join(cv_preds,best_tune_results %>% select(Accession,Observed) %>% distinct())

# cv_preds %<>%
#   mutate(outcome = case_when(Predicted == Observed ~ "Correct",
#                              .default = "Incorrect"))
# 
# # view the incorrect preds for each feature set
# view(cv_preds %>% 
#        filter(feat == "gconly" & outcome == "Incorrect"))
# 
# view(cv_preds %>% 
#        filter(feat == "dintonly" & outcome == "Incorrect"))
# 
# view(cv_preds %>% 
#        filter(feat == "kmer3only" & outcome == "Incorrect"))
# 
# view(cv_preds %>% 
#        filter(feat == "aaconly" & outcome == "Incorrect"))
# 
# view(cv_preds %>% 
#        filter(feat == "dponly" & outcome == "Incorrect"))
# 
# view(cv_preds %>% 
#        filter(feat == "rand" & outcome == "Incorrect"))

best_tune_results %<>%
  mutate(outcome = case_when(Predicted == Observed ~ "Correct",
                             .default = "Incorrect"))

# view the incorrect preds for each feature set
view(best_tune_results %>%
       filter(feat == "gconly" & outcome == "Incorrect"))

view(best_tune_results %>%
       filter(feat == "dintonly" & outcome == "Incorrect"))

view(best_tune_results %>%
       filter(feat == "kmer3only" & outcome == "Incorrect"))

view(best_tune_results %>%
       filter(feat == "aaconly" & outcome == "Incorrect"))

view(best_tune_results %>%
       filter(feat == "dponly" & outcome == "Incorrect"))

view(best_tune_results %>%
       filter(feat == "rand" & outcome == "Incorrect"))