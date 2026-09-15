library(caret)
library(tidyverse)
library(readxl)
library(magrittr)
library(MLeval)
library(writexl)

source("./calc_metrics_GNN_results_loo_funcs.R")

# read in tbl containing subgraph membership information
subgraph_membership <- read_xlsx("../../../../5_phylogeny/3_generate_graph/subgraph_membership_final.xlsx")

# analyse results ---------------------------------------------------------
dirlist <- paste0("./feature_sets/",list.files(path = "./feature_sets",pattern="^results_*"))

cmapply(write_res,results_dir = dirlist)
