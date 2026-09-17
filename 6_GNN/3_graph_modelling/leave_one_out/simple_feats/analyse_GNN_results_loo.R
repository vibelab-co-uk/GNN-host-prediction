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

best_tune <- read_csv("../../kfold_cv_early_stopping/simple_feats/results_summary/best_tune_acc.csv")

# edit best tune tibble to include all rand sets
rand_best_tune <- tibble(nhidden = best_tune$nhidden[best_tune$feat == "rand"],
                         nlayers = best_tune$nlayers[best_tune$feat == "rand"],
                         feat = c("rand1","rand2","rand3","rand4","rand5"),
                         mean_acc = best_tune$mean_acc[best_tune$feat == "rand"],
                         epochs = best_tune$epochs[best_tune$feat == "rand"])

best_tune %<>%
  filter(feat != "rand")

best_tune <- rbind(best_tune,rand_best_tune)

read_best_tunes <- function(feat,nhidden,nlayers){
  blocked_res <- read_csv(paste0("./feature_sets/results_",feat,"/loo_results_",nhidden,"_",nlayers,".csv"))
  blocked_res$feat <- feat
  return(blocked_res)
}

best_tune_results <- do.call(rbind,mapply(read_best_tunes,
                                          feat = best_tune$feat,
                                          nhidden = best_tune$nhidden,
                                          nlayers = best_tune$nlayers,
                                          SIMPLIFY = F))

best_tune_results %<>% 
  mutate(Correct = case_when(Observed == Predicted ~ T,
                             .default = F))

train_incorrect <- best_tune_results %>% 
  filter(Train == T & Correct == F)

nrow(train_incorrect)/nrow(best_tune_results)

# 3.5% of nodes in the train set are incorrectly predicted
pred_node <- best_tune_results %>% 
  filter(Train == F)

view(pred_node %>% 
  count(feat,Correct))

view(pred_node %>% 
  filter(Correct == F) %>% 
  select(Observed,Predicted,prediction_node,feat))
         
         
# read in graph 
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

# hosts
host_labels <- read_csv("../../../../3_host_assignment/all_cluster_hosts.csv")

# read in layout
layout_tbl <- read_csv("../../../../5_phylogeny/3_generate_graph/graph_phylo_plots/final_manual_layout.csv")

layout <- layout_tbl %>% select(x_coord,y_coord)
layout <- as.matrix(layout)

node_id <- inner_join(node_names,host_labels,join_by("Accession" == "ref_accessions"))
V(g)$names <- node_id$Node

# palette <- c(carto_pal(12, "Safe"),"grey1")
palette <- carto_pal(12, "Safe")

mycol = palette[factor(node_id$Host_rank)]

png(filename="../../../../5_phylogeny/3_generate_graph/graph_phylo_plots/host_cols_node_numbers.png",width=1200,heigh=2400,units="px")
plot(g, vertex.size=3,vertex.color=mycol,vertex.label.color="black",vertex.label.degree=pi,vertex.label.dist = 0.5,edge.width=3,layout = layout)
legend('topright',legend=levels(factor((node_id$Host_rank))),col = palette,pch = 16,cex=3)
dev.off()
