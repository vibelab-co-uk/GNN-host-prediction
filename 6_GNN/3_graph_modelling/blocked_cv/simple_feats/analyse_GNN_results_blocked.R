library(tidyverse)
library(readxl)
library(viridis)
library(magrittr)
library(igraph)

`%nin%` <- Negate(`%in%`)

subgraph_id <- read_xlsx("../../../../5_phylogeny/3_generate_graph/subgraph_virus_classes.xlsx")

edgelist <- read_delim("../../../../5_phylogeny/3_generate_graph/graph_data/patristic_paraymxo_graph.edgelist",col_names = F)
edgelist <- edgelist+1

g <- graph_from_edgelist(as.matrix(edgelist),directed = F)

# read in metadata
branch_weights <- read_csv("../../../../5_phylogeny/3_generate_graph/graph_data/patristic_paramyxo_graph_branch_weights.csv")
colnames(branch_weights) <- c("Node","Weight")

node_names <- read_csv("../../../../5_phylogeny/3_generate_graph/graph_data/patristic_paramyxo_graph_node_names.csv")
colnames(node_names) <- c("Node","Accession")

E(g)$weight <- branch_weights$Weight
V(g)$names <- node_names$Accession

# hosts
host_labels <- read_csv("../../../../3_host_assignment/all_cluster_hosts.csv")

# read in layout
layout_tbl <- read_csv("../../../../5_phylogeny/3_generate_graph/graph_phylo_plots/final_manual_layout.csv")

layout <- layout_tbl %>% select(x_coord,y_coord)
layout <- as.matrix(layout)

read_res <- function(results_dir){
  
  res_filename <- list.files(path = results_dir,pattern="^abl_out.*")
  
  # read in full results table
  result_preds <- read_csv(paste0(results_dir,"/",res_filename))
  
  result_preds$feat <- str_split(str_split(results_dir,pattern = "/")[[1]][3],pattern = "_")[[1]][2]
  return(result_preds)
}

dirlist <- paste0("./feature_sets/",list.files(path = "./feature_sets",pattern="^results_*"))

result_preds <- do.call(rbind,lapply(dirlist,read_res))

# process random pseudofeatures
rand_results <- result_preds %>%
  filter(feat %in% c("rand1","rand2","rand3","rand4","rand5"))

# take mean results across the five random pseudofeatures
rand_results %<>%
  group_by(Accession,cv_fold) %>%
  summarise(across(starts_with("tax"),~mean(.x)))

rand_results$feat <- "rand"

feat_results <- result_preds %>%
  filter(feat %nin% c("rand1","rand2","rand3","rand4","rand5")) %>% 
  select(colnames(rand_results))

all_preds <- rbind(feat_results,rand_results)

# make final predictions
all_preds %<>% 
  group_by(Accession,cv_fold,feat) %>% 
  mutate(Prediction_score = max(across(starts_with("tax")))) %>% 
  # subsequent mutate assigns the name of the tax column with the max score (thereby giving the prediction)
  mutate(Predicted = str_sub(colnames(all_preds)[grepl("^tax",colnames(all_preds))][which.max(across(starts_with("tax")))],start = 5))

all_preds$feat <- factor(all_preds$feat,levels = c("gconly","dintonly","kmer3only","aaconly","dponly","rand"))

# bind back in some metadata
all_preds <- left_join(all_preds,result_preds %>% select(Accession,Observed) %>% distinct())

all_preds %<>%
  mutate(outcome = case_when(Predicted == Observed ~ "Correct",
                             .default = "Incorrect"))

palette <- c("wheat3","grey3")

shapes <- c("square","circle")

# we'll now write a function that selects each subgraph and plots
plot_blocked_cv <- function(target_feat,subgraph){
  
  feat_res <- all_preds %>%
    filter(feat == target_feat & cv_fold == subgraph)
  
  train_set <- result_preds %>%
    filter(cv_fold == subgraph & feat == "aaconly") %>% # training set is the same for each cv fold regardless of feat 
    select(Train)
  
  # use a join to order the prediction results correctly
  node_res <- left_join(node_names,feat_res)
  
  V(g)$train <- train_set$Train
  
  V(g)$outcome <- node_res$outcome
  
  mycol <- palette[as.factor(V(g)$outcome)]
  
  myshape <- shapes[as.factor(V(g)$train)]
  
  png(filename=paste0("./results_summary/",target_feat,"/best_tune_",target_feat,"_blocked_cv_",subgraph,"_graph_wlegend.png"),width=2000,heigh=4000,units="px")
  plot(g, vertex.size=3, vertex.label=NA,vertex.color=mycol,vertex.shape=myshape,layout = layout,edge.width=3)
  legend('topright',legend=levels(factor((V(g)$outcome))),col = palette,pch = 16,cex=3)
  legend('bottomright',legend = c("Train","Test"),pch = c(1,0),cex=3)
  dev.off()
}

# define a function to apply function to every combination of 2 arguments
cmapply <- function(FUN, ..., MoreArgs = NULL){  
  # expand a grid of all argument combinations  
  l <- expand.grid(..., stringsAsFactors=FALSE)    
  # apply the function  
  .mapply(FUN=FUN, dots=unname(l), MoreArgs = MoreArgs)}

cmapply(plot_blocked_cv,target_feat = unique(all_preds$feat),subgraph = unique(all_preds$cv_fold))
