library(tidyverse)
library(readxl)
library(viridis)
library(magrittr)
library(igraph)
library(rcartocolor)
library(RColorBrewer)
library(terra)
library(ape)
library(ggtree)


# define a function to apply function to every combination of 2 arguments
cmapply <- function(FUN, ..., MoreArgs = NULL){  
  # expand a grid of all argument combinations  
  l <- expand.grid(..., stringsAsFactors=FALSE)    
  # apply the function  
  .mapply(FUN=FUN, dots=unname(l), MoreArgs = MoreArgs)}

best_tune <- read_csv("../../kfold_cv_early_stopping/simple_feats/results_summary/best_tune_acc.csv")

rand_best_tune <- tibble(nhidden = best_tune$nhidden[best_tune$feat == "rand"],
                         nlayers = best_tune$nlayers[best_tune$feat == "rand"],
                         feat = c("rand1","rand2","rand3","rand4","rand5"),
                         mean_acc = best_tune$mean_acc[best_tune$feat == "rand"],
                         epochs = best_tune$epochs[best_tune$feat == "rand"])

best_tune %<>%
  filter(feat != "rand")

best_tune <- rbind(best_tune,rand_best_tune)

read_best_tunes <- function(feat,nhidden,nlayers){
  blocked_res <- read_csv(paste0("./feature_sets/results_",feat,"/loo_results_",nhidden,"_",nlayers,"convlayers_weighted.csv"))
  blocked_res$feat <- feat

  return(blocked_res)
}

result_preds <- do.call(rbind,mapply(read_best_tunes,
                                          feat = best_tune$feat,
                                          nhidden = best_tune$nhidden,
                                          nlayers = best_tune$nlayers,
                                          SIMPLIFY = F,
                                          USE.NAMES = F))

result_preds %>% 
  count(feat,Correct)
         
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

# plot with nodes coloured by prediction score
# make a vector of prediction score for the observed host
result_preds$observed_colnames <- paste0("tax_",result_preds$Observed)

result_preds %<>%
  rowwise %>% 
  mutate(true_p_score = cur_data()[[observed_colnames]]) %>% 
  select(-observed_colnames)

plot_feat_pscore_graph <- function(target_feat){
  
  palette_order_correct <- result_preds %>% 
    filter(feat == target_feat & Correct == T) %>% 
    arrange(true_p_score) %>% 
    select(true_p_score,Accession)
  
  correct_pal <- colorRampPalette(c("white","steelblue2","steelblue3"))(round(nrow(palette_order_correct)*1.2))
  
  correct_pal <- correct_pal[(round(nrow(palette_order_correct)*0.2)+1):length(correct_pal)]
  
  palette_order_correct$col_pal <- correct_pal
  
  palette_order_incorrect <- result_preds %>% 
    filter(feat == target_feat & Correct == F) %>% 
    arrange(-true_p_score) %>% 
    select(true_p_score,Accession)
  
  incorrect_pal <- colorRampPalette(c("white","firebrick2","firebrick"))(round(nrow(palette_order_incorrect)*1.2))
  
  incorrect_pal <- incorrect_pal[(round(nrow(palette_order_incorrect)*0.2)+1):length(incorrect_pal)]
  
  palette_order_incorrect$col_pal <- incorrect_pal
  
  palette_order <- rbind(palette_order_correct,palette_order_incorrect)
  
  node_id_cols <- left_join(node_id,palette_order)
  
  hinge <- min(result_preds %>% 
                 filter(feat == target_feat & Correct == T) %>% 
                 select(true_p_score))
  
  hinge <- hinge*100
  
  correct_legend_pal <- colorRampPalette(c("white","steelblue2","steelblue","steelblue4"))((100 - round(hinge))*1.2)
  
  correct_legend_pal <- correct_legend_pal[round(((100 - round(hinge))*0.2)+1):length(correct_legend_pal)]
  
  incorrect_legend_pal <- colorRampPalette(c("white","firebrick2","firebrick","firebrick4"))(round(hinge)*1.2)
  
  incorrect_legend_pal <- incorrect_legend_pal[round((round(hinge)*0.2)+1):length(incorrect_legend_pal)]
  
  pdf(file=paste0("./results_summary/prediction_cols_graph_plot_",target_feat,".pdf"),width=24,height=48)
  plot(g, vertex.size=3,vertex.color=node_id_cols$col_pal,vertex.label=NA,edge.width=4,layout = layout)
  legend_cont(x=-0.7,y=-1.4,horiz=T,legend=c(0,1),col = c(rev(incorrect_legend_pal),correct_legend_pal), size = c(1,6),
              title = "Correct label probability",at = c(0,0.25,0.5,0.75,1),cex=4,title.cex=5)
  dev.off()
}

lapply(unique(result_preds$feat),plot_feat_pscore_graph)

# generate a phylogeny showing predictions across all labels
# read in BEAST phylogeny
p_phylo <- ape::read.nexus("../../../../5_phylogeny/2_beast_phylo/GTR_HMC_HPSTR_tree_50.nxs")

hosts <- read_csv("../../../../3_host_assignment/all_cluster_hosts.csv")

# create tibble of tip labels and bind in hosts
phylo_tips <- tibble(ref_accessions = p_phylo$tip.label)

tip_hosts <- left_join(phylo_tips, hosts)

tip_hosts %<>%
  mutate(Host_rank = case_when(is.na(Host_rank) == T ~ "Undetermined",
                               .default = Host_rank))

# label the phylogeny with hosts and drop tips with undetermined hosts
p_phylo$tip.label <- factor(tip_hosts$Host_rank,levels=c("Primates","Rodentia","Carnivora","Chiroptera","Artiodactyla","Perissodactyla","Scandentia","Pholidota","Eulipotyphla",
                                                         "Aves","Reptilia","Fish","Undetermined"))

p_phylo <- ape::drop.tip(p_phylo,"Undetermined")

tip_hosts %<>%
  filter(Host_rank != "Undetermined")


# set up a colour palette for plotting
palette <- c(carto_pal(12, "Safe"),"grey1")
mycol = palette[factor(tip_hosts$Host_rank,levels=c("Primates","Rodentia","Carnivora","Chiroptera","Artiodactyla","Perissodactyla","Scandentia","Pholidota","Eulipotyphla",
                                                    "Aves","Reptilia","Fish"))]

# label the phylogeny with hosts
p_phylo$tip.label <- tip_hosts$ref_accessions

p_phylo <- left_join(p_phylo,tip_hosts,join_by("label" == "ref_accessions"))

p_point <- ggtree(p_phylo,ladderize=F) +
  geom_tippoint(size = 1.5,color=mycol)

p_accession <- ggtree(p_phylo,ladderize=F) +
  geom_tiplab(size = 1.5,color=mycol)

plot_feat_pscore_phylo <- function(target_feat){

  plot_res <- result_preds %>% 
    filter(feat == target_feat) %>% 
    select(Accession,starts_with("tax_"))
  
  plot_res <- data.frame(plot_res)
  
  rownames(plot_res) <- plot_res$Accession
  
  plot_res <- plot_res[,-1]
  
  gheatmap(p_point, plot_res, offset=0.1, width=0.4,colnames_angle=-45,hjust = 0,font.size = 3,
           custom_column_labels=str_sub(colnames(plot_res),start=5)) +
    scale_fill_continuous(palette = c("firebrick4","firebrick3","firebrick2","steelblue2","steelblue","steelblue4"),name = "Output probability") +
    ggtree::vexpand(.2, -1) +
    coord_cartesian(clip = "off")
  
  ggsave("new_seqs_heatmap_tree.pdf",width = 5000,height = 8000,units = "px",dpi=600)
  
} 

plot_res2 <- result_preds %>% 
  filter(feat == target_feat) %>% 
  select(Accession,true_p_score)

plot_res2 <- data.frame(plot_res2)

res_rownames <- plot_res2$Accession

plot_res2 <- plot_res2[,-1]

plot_res2 <- data.frame(plot_res2)

rownames(plot_res2) <- res_rownames

gheatmap(p_point, plot_res2, offset=0.1, width=0.4,colnames_angle=-45,hjust = 0,font.size = 3,
         custom_column_labels=str_sub(colnames(plot_res),start=5)) +
  scale_fill_continuous(palette = c("white","steelblue2","steelblue","steelblue4"),name = "Output probability") +
  ggtree::vexpand(.2, -1) +
  coord_cartesian(clip = "off")
