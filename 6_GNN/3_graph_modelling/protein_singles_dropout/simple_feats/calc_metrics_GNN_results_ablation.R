library(caret)
library(tidyverse)
library(readxl)
library(magrittr)
library(ranger)
library(MLeval)
library(writexl)

source("./calc_metrics_GNN_results_funcs.R")

# read in tbl containing subgraph membership information
subgraph_membership <- read_xlsx("../../../../5_phylogeny/3_generate_graph/subgraph_membership_final.xlsx")

# analyse results ---------------------------------------------------------
dirlist <- paste0("./feature_sets/",list.files(path = "./feature_sets",pattern="^results_*"))

lapply(dirlist,write_res)

# assess performance over different feature sets --------------------------
all_agg <- do.call(rbind,lapply(dirlist,cat_agg,hidden = nhidden,model = model_params))
all_roc <- do.call(rbind,lapply(dirlist,cat_roc,hidden = nhidden,model = model_params))
all_noinfo <- do.call(rbind,lapply(dirlist,cat_noinfo,hidden = nhidden,model = model_params))

all_agg$Ablation <- as.factor(all_agg$Ablation)

subgraph_agg <- all_agg %>% 
  group_by(Ablation,Subgraph,feat,tax_level) %>% 
  dplyr::summarise(mean_accuracy = mean(Accuracy),mean_macro = mean(macro_f1,na.rm=T))

all_roc %<>% 
  pivot_wider(names_from = Metric,values_from = c(Score,CI))

subgraph_roc <- all_roc %>% 
  group_by(Ablation,feat,Subgraph) %>% 
  dplyr::summarise(mean_ROC = mean(`Score_AUC-ROC`,na.rm=T),mean_PRG = mean(`Score_AUC-PRG`,na.rm=T))

subgraph_noinfo <- all_noinfo %>% 
  group_by(Ablation,Subgraph,feat,tax_level) %>% 
  dplyr::summarise(mean_accuracy = mean(Accuracy),mean_macro = mean(macro_f1,na.rm=T))

# write outputs 
write_xlsx(subgraph_agg, paste0("./results_summary/subgraph_agg_performance.xlsx"))
write_xlsx(subgraph_roc, paste0("./results_summary/subgraph_roc_performance.xlsx"))
write_xlsx(subgraph_noinfo, paste0("./results_summary/subgraph_noinfo_performance.xlsx"))

write_xlsx(all_agg, paste0("./results_summary/agg_performance_full.xlsx"))
write_xlsx(all_roc, paste0("./results_summary/roc_performance_full.xlsx"))
write_xlsx(all_noinfo, paste0("./results_summary/noinfo_performance_full.xlsx"))



