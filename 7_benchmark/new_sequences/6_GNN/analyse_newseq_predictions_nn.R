library(tidyverse)
library(magrittr)
library(ape)
library(ggtree)
library(rcartocolor)

# read in index of new sequence IDs
new_seqs_hosts <- read_csv("../3_host_assignment/all_cluster_hosts.csv")
new_seqs_hosts <- as.vector(new_seqs_hosts$ref_accessions)

read_res <- function(x){
  
  results <- read_csv(paste0(x,"/new_seqs_output.csv"))
  
  results$feat <- str_split_i(x,pattern = "_",i = 2)
  
  results %<>%
    mutate(new_seq = case_when(Accession %in% new_seqs_hosts ~ T,
                               .default = F)) %>% 
    mutate(Outcome = case_when(Observed == Predicted ~ "Correct",
                               .default = "Incorrect"))
  
  return(results)
}

dirlist_gnn <- list.files(path = ".",pattern = "results_*")

gnn_results <- do.call(rbind,lapply(dirlist_gnn,read_res))

gnn_results %<>% 
  filter(feat == "dponly")

gnn_results$feat <- "dponly_gnn"

# quick look at results
gnn_summary <- gnn_results %>% 
  filter(new_seq == T)

gnn_summary %>% 
  count(feat,Outcome)
# read in results from non-graph NN
dirlist_nn <- list.files(path = "../../nn_benchmark/dense_nn/new_sequences",pattern = "results_*")

nn_results <- do.call(rbind,lapply(dirlist_nn,read_res))

nn_results %<>% 
  filter(feat == "dponly")

nn_results$feat <- "dponly_nn"

# quick look at results
nn_summary <- nn_results %>% 
  filter(new_seq == T)

nn_summary %>% 
  count(feat,Outcome)

all_results <- rbind(gnn_results,nn_results)

# read in tree with all sequences
# read in maximum likelihood phylogeny
p_phylo <- ape::read.tree("../../new_sequences/5_phylogeny/2_iqtree_phylogeny/new_seqs.treefile")
old_hosts <- read_csv("../../../3_host_assignment/all_cluster_hosts.csv")
new_hosts <- read_csv("../../new_sequences/3_host_assignment/all_cluster_hosts.csv")

all_hosts <- rbind(old_hosts,new_hosts)

# root tree using the highly divergent seal-associated virus
p_phylo <- root(p_phylo, outgroup = "PX841364.1",resolve.root = T)

# create tibble of tip labels and bind in hosts
phylo_tips <- tibble(ref_accessions = p_phylo$tip.label)

tip_hosts <- left_join(phylo_tips, all_hosts)

tip_hosts %<>%
  drop_na()

# set up a colour palette for plotting
palette <- c(carto_pal(12, "Safe"),"grey1")
mycol = palette[factor(tip_hosts$Host_rank,levels=c("Primates","Rodentia","Carnivora","Chiroptera","Artiodactyla","Perissodactyla","Scandentia","Pholidota","Eulipotyphla",
                                                    "Aves","Reptilia","Fish"))]
# label the phylogeny with hosts
p_phylo$tip.label <- tip_hosts$ref_accessions

# axis_fun <- function() {
#   col <- "black"
#   for (i in 1:2)
#     axis(i, col = col, col.ticks = col, col.axis = col, las = 1)
#   box(lty = "19")
# }
# 
# plot(p_phylo, tip.color = mycol); axis_fun()

# read in BLAST results
blast_res_new <- read_csv("../../blastn_benchmark/new_sequences/new_seqs_blast_pred_results.csv")
blast_res_old <- read_csv("../../blastn_benchmark/whole_graph/blast_pred_results.csv")

blast_res <- rbind(blast_res_new,blast_res_old)

blast_res$feat <- "BLAST"

blast_res %<>%
  select(Accession,feat,Outcome)

# drop blast results for unknown hosts
blast_res %<>%
  filter(Accession %in% all_results$Accession) %>% 
  distinct()

plot_res <- all_results %>% 
  select(Accession,feat,Outcome)

plot_res <- rbind(plot_res,blast_res)

plot_res %<>%
  mutate(Outcome = case_when(Accession %in% new_hosts$ref_accessions ~ Outcome,
                             .default = "Train"))

plot_res %<>%
  pivot_wider(names_from = feat,values_from = Outcome,names_sep = " ")

plot_res <- data.frame(plot_res)

rownames(plot_res) <- paste0(plot_res$Accession,plot_res$model)

plot_res <- plot_res[,-1]

p_phylo <- left_join(p_phylo,tip_hosts,join_by("label" == "ref_accessions"))

p_point <- ggtree(p_phylo,ladderize=F) +
  geom_tippoint(size = 1.5,color=mycol)

p_accession <- ggtree(p_phylo,ladderize=F) +
  geom_tiplab(size = 1.5,color=mycol)

gheatmap(p_point, plot_res, offset=0.1, width=0.4, legend_title="genotype",colnames_angle=-45,hjust = 0,colnames_offset_y = 5,font.size = 3,
         custom_column_labels=c("Dipeptide composition - GNN",
                                "Dipeptide composition - NN",
                                "BLAST benchmark")) +
  scale_fill_manual(values = c("steelblue2","firebrick2","white"),labels = c("Correct","Incorrect","Train"),name = NULL) +
  ggtree::vexpand(.2, -1) +
  coord_cartesian(clip = "off")

ggsave("new_seqs_heatmap_tree.png",width = 5000,height = 8000,units = "px",dpi=600)

gheatmap(p_accession, plot_res, offset=0.4, width=0.4, legend_title="genotype",colnames_angle=-45,hjust = 0,colnames_offset_y = 5,font.size = 3,
         custom_column_labels=c("Dipeptide composition - GNN",
                                "Dipeptide composition - NN",
                                "BLAST benchmark")) +
  scale_fill_manual(values = c("steelblue2","firebrick2","white"),labels = c("Correct","Incorrect","Train"),name = NULL) +
  ggtree::vexpand(.2, -1) +
  coord_cartesian(clip = "off")

ggsave("new_seqs_accessions_heatmap_tree.png",width = 5000,height = 8000,units = "px",dpi=600)

# reduced version
# add dummy columns for heatmap spacing
plot_res$dummy1 <- "Train"
plot_res$dummy2 <- "Train"

short_res <- plot_res %>% 
  select(dponly_gnn,dummy1,dponly_nn,dummy2,BLAST)

gheatmap(p_point, short_res, offset=0.1, width=0.3, legend_title="genotype",colnames_angle=-25,hjust = 0,colnames_offset_y = 5,font.size = 3,
         custom_column_labels=c("Dipeptide composition - GNN",
                                "",
                                "Dipeptide composition - NN",
                                "",
                                "BLAST benchmark")) +
  scale_fill_manual(values = c("steelblue2","firebrick2","white"),labels = c("Correct","Incorrect","Train"),name = NULL) +
  ggtree::vexpand(.2, -1) +
  coord_cartesian(clip = "off")

ggsave("new_seqs_short_heatmap_tree_nn.svg",width = 5000,height = 8000,units = "px",dpi=600)

