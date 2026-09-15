library(tidyverse)
library(readxl)
library(writexl)
library(caret)

### this script generates train/test sets for a range of validation methods
# in all sample sets, numbered indices represent nodes for use in TEST sets

### define functions
`%nin%` <- Negate(`%in%`)

# define a function to apply function to every combination of arguments
cmapply <- function(FUN, ..., MoreArgs = NULL){  
  
  l <- expand.grid(..., stringsAsFactors=FALSE) # expand a grid of all argument combinations
  
  .mapply(FUN=FUN, dots=unname(l), MoreArgs = MoreArgs)}   # apply the function  

# read in list of graph nodes and metadata
nodes <- read_csv("../../5_phylogeny/3_generate_graph/graph_data/patristic_paramyxo_graph_node_names.csv")
colnames(nodes) <- c("index","Accession")

# host labels
host_labels <- read_csv("../../3_host_assignment/all_cluster_hosts.csv")

nodes <- left_join(nodes,host_labels,join_by("Accession" == "ref_accessions"))

# convert to factor so createFolds function will attempt to balance class distribution between within splits 
nodes$Host_rank <- factor(nodes$Host_rank)

# subgraph membership
subgraphs <- read_xlsx("../../5_phylogeny/3_generate_graph/subgraph_membership_final.xlsx")

nodes <- left_join(nodes,subgraphs)

# generate kfold cross validation samples ----------------------------------
cv_k <- 5 # set desired number of kfolds

# create cross validation folds
cvFolds <- createFolds(nodes$Host_rank, k=cv_k)

# pad all folds to the same length
max_l <- max(c(length(cvFolds$Fold1),length(cvFolds$Fold2),length(cvFolds$Fold3),length(cvFolds$Fold4),length(cvFolds$Fold5)))

length(cvFolds$Fold1) <- max_l
length(cvFolds$Fold2) <- max_l
length(cvFolds$Fold3) <- max_l
length(cvFolds$Fold4) <- max_l
length(cvFolds$Fold5) <- max_l

# bind into a dataframe
folds_df <- as.data.frame(do.call(cbind,cvFolds))

# write outputs
write_csv(folds_df,paste0("cv_folds_",cv_k,".csv"))

# generate another set of kfolds for initial hyperparameter tuning
# create cross validation folds
cvFolds <- createFolds(nodes$Host_rank, k=cv_k)

# pad all folds to the same length
max_l <- max(c(length(cvFolds$Fold1),length(cvFolds$Fold2),length(cvFolds$Fold3),length(cvFolds$Fold4),length(cvFolds$Fold5)))

length(cvFolds$Fold1) <- max_l
length(cvFolds$Fold2) <- max_l
length(cvFolds$Fold3) <- max_l
length(cvFolds$Fold4) <- max_l
length(cvFolds$Fold5) <- max_l

# bind into a dataframe
folds_df <- as.data.frame(do.call(cbind,cvFolds))

# write outputs
write_csv(folds_df,paste0("cv_folds_",cv_k,"_tuning.csv"))

# generate blocked cross validation samples-------------------------------
# write a vector of subgraphs we want to include (very small subgraphs with few sequences are excluded)
target_subgraphs <- c(1,2,3,4,6,8,9,10)

blocked_folds_df <- data.frame(do.call(cbind,rep(nodes %>% select(index),each=length(target_subgraphs))))

colnames(blocked_folds_df) <- paste0("subgraph",target_subgraphs)

# set nodes not in the current subgraph to NA
for(i in 1:length(target_subgraphs)){
  blocked_folds_df[,i][nodes$Subgraph != str_extract(colnames(blocked_folds_df)[i],"\\d+")] <- NA # here we use the column name to filter for nodes in the correct subgraph
}

write_csv(blocked_folds_df,paste0("blocked_cv_folds.csv"))

# generate ensemble folds outside of each block for ensemble modelling
# repeat inverse of previous operation to give numbered indices of nodes not in each subgraph
ensemble_folds_df <- data.frame(do.call(cbind,rep(nodes %>% select(index),each=length(target_subgraphs))))

colnames(ensemble_folds_df) <- paste0("subgraph",target_subgraphs)

# set nodes not in the current subgraph to NA
for(i in 1:length(target_subgraphs)){
  ensemble_folds_df[,i][nodes$Subgraph == str_extract(colnames(ensemble_folds_df)[i],"\\d+")] <- NA # here we use the column name to filter for nodes in the correct subgraph
}

# for each subgraph block, generate 5 ensemble folds from outside the subgraph
for(i in 1:ncol(ensemble_folds_df)){
  ensemble_folds <- createFolds(nodes$Host_rank[ensemble_folds_df[,1]],k = cv_k) # create folds will create an additional fold containing the NAs (Fold1), we can drop this
  
  # pad all folds to the same length
  max_l <- max(c(length(ensemble_folds$Fold2),length(ensemble_folds$Fold3),length(ensemble_folds$Fold4),length(ensemble_folds$Fold5),length(ensemble_folds$Fold6)))
  
  length(ensemble_folds$Fold2) <- max_l
  length(ensemble_folds$Fold3) <- max_l
  length(ensemble_folds$Fold4) <- max_l
  length(ensemble_folds$Fold5) <- max_l
  length(ensemble_folds$Fold6) <- max_l
  
  # bind into a dataframe
  folds_df <- as.data.frame(do.call(cbind,ensemble_folds[-1])) # bind together the folds not containing the NAs
  colnames(folds_df) <- c("Fold1","Fold2","Fold3","Fold4","Fold5")
  
  # write outputs
  write_csv(folds_df,paste0("ensemble_folds_subgraph",str_extract(colnames(blocked_folds_df)[i],"\\d+"),".csv"))
}

# generate blocked cross validation samples with reintroductions
reintro_blocked_folds_df <- blocked_folds_df

sample_subgraphs <- function(n,iteration){
  
  set.seed(iteration*123) # set a random seed based on the iteration
  
  # randomly sample from within the target subgraph
  for(i in 1:length(target_subgraphs)){
  sub_target <- nodes %>% 
    filter(Subgraph == str_extract(colnames(reintro_blocked_folds_df)[i],"\\d+")) %>% 
    slice_sample(n=n)
  
  reintro_blocked_folds_df[,i][sub_target$index] <- NA # replace indices of reintroduced sequences with NAs (these will no longer be in the test set)
  }
  
  write_csv(reintro_blocked_folds_df,paste0("reintro_blocked_cv_folds_n",n,"_iteration",iteration,".csv"))
}

# set the number of random reintroductions per subgraph and the number of repeats of random selection
cmapply(sample_subgraphs,n = 1:5, iteration = 1:5) 
