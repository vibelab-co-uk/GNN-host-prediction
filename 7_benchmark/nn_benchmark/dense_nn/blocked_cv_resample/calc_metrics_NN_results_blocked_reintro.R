library(caret)
library(tidyverse)
library(readxl)
library(magrittr)
library(ranger)
library(MLeval)
library(writexl)

source("./calc_metrics_NN_results_reintro_funcs.R")

reintro = c(1,2,3,4,5)
iter = c(1,2,3,4,5)

# read in tbl containing subgraph membership information
subgraph_membership <- read_xlsx("../../../../5_phylogeny/3_generate_graph/subgraph_membership_final.xlsx")

# analyse results ---------------------------------------------------------
dirlist <- paste0("./feature_sets/",list.files(path = "./feature_sets",pattern="^results_*"))

cmapply(write_res,results_dir = dirlist,reintro = reintro,iter = iter)

# assess performance over different feature sets --------------------------
all_agg <- do.call(rbind,cmapply(cat_agg,results_dir = dirlist,reintro = reintro,iter = iter))
all_roc <- do.call(rbind,cmapply(cat_roc,results_dir = dirlist,reintro = reintro,iter = iter))
all_noinfo <- do.call(rbind,cmapply(cat_noinfo,results_dir = dirlist,reintro = reintro,iter = iter))

subgraph_agg <- all_agg %>% 
  group_by(Block,feat,tax_level,N_samples) %>% 
  dplyr::summarise(mean_accuracy = mean(Accuracy),mean_macro = mean(macro_f1,na.rm=T),sd_accuracy = sd(Accuracy))

all_roc %<>% 
  pivot_wider(names_from = Metric,values_from = c(Score,CI))

subgraph_roc <- all_roc %>% 
  group_by(Block,feat,N_samples) %>% 
  dplyr::summarise(mean_ROC = mean(`Score_AUC-ROC`,na.rm=T),mean_PRG = mean(`Score_AUC-PRG`,na.rm=T))

subgraph_noinfo <- all_noinfo %>% 
  group_by(Block,feat,tax_level,N_samples) %>% 
  dplyr::summarise(mean_accuracy = mean(Accuracy),mean_macro = mean(macro_f1,na.rm=T))

# write outputs 
setwd("./results_summary")

write_xlsx(subgraph_agg, paste0("subgraph_agg_performance_best_tune.xlsx"))
write_xlsx(subgraph_roc, paste0("subgraph_roc_performance_best_tune.xlsx"))
write_xlsx(subgraph_noinfo, paste0("subgraph_noinfo_performance_best_tune.xlsx"))

write_xlsx(all_agg, paste0("agg_performance_full_best_tune.xlsx"))
write_xlsx(all_roc, paste0("roc_performance_full_best_tune.xlsx"))
write_xlsx(all_noinfo, paste0("noinfo_performance_full_best_tune.xlsx"))



