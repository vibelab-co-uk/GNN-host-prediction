library(caret)
library(tidyverse)
library(readxl)
library(magrittr)
library(MLeval)
library(writexl)
library(RColorBrewer)
library(philentropy)
library(vegan)
library(rcartocolor)
library(igraph)

source("./embedding_analysis_funcs.R")

`%nin%` <- Negate(`%in%`)

epochs <- seq(from=0,to=500,by=25) # epochs for tsne plots
cv_folds <- c(1,2,3,4,5)

dirlist <- list.files(path = "./feature_sets",pattern="^results_*")

# combine result dumpouts into single file ----------------------------------
host_labels <- lapply(dirlist,combine_res)

host_labels <- host_labels[[1]] # labels are identical for all feature sets

#read in the phylogeny 
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

# relabel tips with accessions instead of hosts
p_phylo$tip.label <- tip_hosts$ref_accessions

# plot tsne ---------------------------------------------------------
tsne_files <- lapply(X = paste0("feature_sets/",dirlist),FUN = list.files, pattern = "tsne_*")

all_tsne <- do.call(rbind,mapply(read_tsne,file = tsne_files, dir = dirlist,SIMPLIFY= F))

colnames(all_tsne) <- c("tSNE1","tSNE2","cv_fold","Accession","Observed","epoch","feat")

pal <- carto_pal(12, "Safe")

all_tsne$epoch <- as.numeric(all_tsne$epoch)

lapply(all_tsne %>% group_by(feat,cv_fold) %>% group_split(),tsne_plots)

### embedding/patristic distance correlation

# calculate patristic distance --------------------------------------------

pat_distances <- ape::cophenetic.phylo(p_phylo)

# normalise 
pat_distances <- (pat_distances-min(pat_distances))/(max(pat_distances)-min(pat_distances))

# process embeddings ------------------------------------------------------
embedding_files <- lapply(X = paste0("feature_sets/",dirlist),FUN = list.files, pattern = "raw_embeddings_*")

all_embedding_dist <- do.call(rbind,mapply(pat_embedding_similarity,file = embedding_files, dir = dirlist,SIMPLIFY= F))

# take average of all metrics across cv folds
all_embedding_dist %<>%
  group_by(feature_set,epoch) %>% 
  summarise(sum_residuals = mean(sum_residuals),sum_square_residuals = mean(sum_square_residuals),L2_norm = mean(L2_norm))

all_embedding_dist$epoch <- as.numeric(all_embedding_dist$epoch)

### process and visualise distance data
write_xlsx(all_embedding_dist,"./results_summary/embeddings_patristic_metrics_table.xlsx")

# all_embedding_dist <- read_xlsx("./results_summary/embeddings_patristic_metrics_table.xlsx")

# read in loss and accuracy information from the training logs ----------
log_files <- lapply(X = paste0("feature_sets/",dirlist),FUN = list.files, pattern = "training_log_*")

all_training_logs <- do.call(rbind,mapply(FUN = read_acc_loss,files = log_files,dir = dirlist,SIMPLIFY= F))

all_training_logs$epoch <- as.numeric(all_training_logs$epoch)
all_training_logs$training_loss <- as.numeric(all_training_logs$training_loss)
all_training_logs$training_acc <- as.numeric(all_training_logs$training_acc)
all_training_logs$test_loss <- as.numeric(all_training_logs$test_loss)
all_training_logs$test_acc <- as.numeric(all_training_logs$test_acc)

# take mean results across CV folds for the training logs
all_training_logs %<>% 
  group_by(feat,epoch) %>% 
  summarise(training_loss = mean(training_loss),training_acc = mean(training_acc),test_loss = mean(test_loss),test_acc = mean(test_acc))

# generate subsets of only the random pseudofeatures
rand_error <- all_embedding_dist %>% 
  filter(feature_set %in% c("rand1","rand2","rand3","rand4","rand5"))

# generate subsets of only the random pseudofeatures
rand_training_logs <- all_training_logs %>% 
  filter(feat %in% c("rand1","rand2","rand3","rand4","rand5"))

### generate normalised features for plotting against accuracy
# normalise model loss to 0-100
all_training_logs$training_loss_norm <- ((all_training_logs$training_loss-min(all_training_logs$training_loss))/(max(all_training_logs$training_loss)-min(all_training_logs$training_loss)))*100
all_training_logs$test_loss_norm <- ((all_training_logs$test_loss-min(all_training_logs$test_loss))/(max(all_training_logs$test_loss)-min(all_training_logs$test_loss)))*100

# normalise distance metrics to 0-100 for plotting against accuracy
all_embedding_dist$sum_residuals_norm <- ((all_embedding_dist$sum_residuals-min(all_embedding_dist$sum_residuals))/(max(all_embedding_dist$sum_residuals)-min(all_embedding_dist$sum_residuals)))*100
all_embedding_dist$sum_square_residuals_norm <- ((all_embedding_dist$sum_square_residuals-min(all_embedding_dist$sum_square_residuals))/(max(all_embedding_dist$sum_square_residuals)-min(all_embedding_dist$sum_square_residuals)))*100
# all_embedding_dist$mantel_correlation_norm <- ((all_embedding_dist$mantel_correlation-min(all_embedding_dist$mantel_correlation))/(max(all_embedding_dist$mantel_correlation)-min(all_embedding_dist$mantel_correlation)))*100

# note that for random features we normalise using the full range of values to ensure consistency with other features
rand_error$sum_square_norm <- ((rand_error$sum_square_residuals-min(all_embedding_dist$sum_square_residuals))/(max(all_embedding_dist$sum_square_residuals)-min(all_embedding_dist$sum_square_residuals)))*100

rand_training_logs$training_loss_norm <- ((rand_training_logs$training_loss-min(all_training_logs$training_loss))/(max(all_training_logs$training_loss)-min(all_training_logs$training_loss)))*100
rand_training_logs$test_loss_norm <- ((rand_training_logs$test_loss-min(all_training_logs$test_loss))/(max(all_training_logs$test_loss)-min(all_training_logs$test_loss)))*100

# take mean results across the five random pseudofeatures
rand_training_logs <- rand_training_logs %>%
  group_by(epoch) %>%
  summarise(training_acc_sd = sd(training_acc),
            training_acc_mean = mean(training_acc),
            test_acc_sd = sd(test_acc),
            test_acc_mean = mean(test_acc),
            training_loss_sd = sd(training_loss),
            training_loss_mean = mean(training_loss),
            test_loss_sd = sd(test_loss),
            test_loss_mean = mean(test_loss))

rand_training_logs$feature_set <- factor("rand",levels = c("dponly","rand"))

rand_error <- rand_error %>%
  group_by(epoch) %>%
  summarise(sum_square_sd = sd(sum_square_residuals),
            sum_square_mean = mean(sum_square_residuals),
            # mantel_signif_mean = mean(mantel_significance),
            # mantel_cor_mean = mean(mantel_correlation),
            sum_square_sd_norm = sd(sum_square_norm),
            sum_square_mean_norm = mean(sum_square_norm))

# take the average of training and test accuracy to give whole graph accuracy
all_training_logs$all_acc_mean <- rowMeans(all_training_logs %>% ungroup() %>%  select(test_acc,training_acc))

rand_training_logs$all_acc_mean <- rowMeans(rand_training_logs %>% ungroup() %>%  select(test_acc_mean,training_acc_mean))

rand_training_logs$all_acc_sd <- rowMeans(rand_training_logs %>% ungroup() %>%  select(test_acc_sd,training_acc_sd))


rand_error$feature_set <- factor("rand",levels = c("dponly","rand"))

all_embedding_dist_norand <- all_embedding_dist %>%
  filter(feature_set %nin% c("rand1","rand2","rand3","rand4","rand5"))

all_training_logs_norand <- all_training_logs %>%
  filter(feat %nin% c("rand1","rand2","rand3","rand4","rand5"))

# process data for plotting
all_embedding_dist_norand$feature_set <- factor(all_embedding_dist_norand$feature_set,levels = c("dponly"))

all_training_logs_norand$feat <- factor(all_training_logs_norand$feat,levels = c("dponly"))


palette_OkabeIto <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                      "#0072B2", "#D55E00", "#CC79A7", "#999999")


palette_OkabeIto_2 <- c("#CC79A7", "#999999")


# plot output
# ggplot(all_embedding_dist,aes(x = epoch,y = sum_residuals,col=feature_set)) +
#   geom_point()+
#   geom_line()+
#   scale_colour_manual(labels = c("GC %","Dinucleotide composition","3mer composition","Amino acid composition",
#                                  "Dipeptide composition","Random pseudofeature 1","Random pseudofeature 2",
#                                  "Random pseudofeature 3","Random pseudofeature 4","Random pseudofeature 5"),values = palette_OkabeIto_6)+
#   theme_bw(base_size = 12)+
#   xlab("Training epoch")+
#   ylab("Sum of residuals")+
#   guides(colour = guide_legend(title = "Feature set",override.aes = list(alpha = 1)))+
#   facet_wrap(~model_params,ncol = 3,labeller = facet_titles) +
#   theme(strip.background = element_rect(fill="white"))
# 
# ggsave("sum_residuals_plot.png",height = 2440,width = 2880,unit="px")

ggplot(all_embedding_dist_norand,aes(x = epoch,y = sum_square_residuals,col=feature_set)) +
  geom_point()+
  geom_line()+
  geom_pointrange(data = rand_error,aes(x = epoch,y=sum_square_mean,ymin = sum_square_mean-sum_square_sd,ymax = sum_square_mean+sum_square_sd),size = 0.2)+
  geom_line(data = rand_error,aes(x = epoch,y=sum_square_mean))+
  scale_colour_manual(labels = c("Dipeptide composition","Random pseudofeature"),values = palette_OkabeIto_2)+
  theme_bw(base_size = 12)+
  xlab("Training epoch")+
  ylab("Sum of squared residuals")+
  guides(colour = guide_legend(title = "Feature set",override.aes = list(alpha = 1)))+
  # facet_wrap(~model_params,ncol = 3,labeller = facet_titles) +
  theme(strip.background = element_rect(fill="white"))

ggsave("./results_summary/sum_square_residuals_plot.png",height = 2440,width = 2880,unit="px")

# ggplot(all_embedding_dist,aes(x = epoch,y = L2_norm,col=feature_set)) +
#   geom_point()+
#   geom_line()+
#   scale_colour_manual(labels = c("GC %","Dinucleotide composition","3mer composition","Amino acid composition",
#                                  "Dipeptide composition","Random pseudofeature 1","Random pseudofeature 2",
#                                  "Random pseudofeature 3","Random pseudofeature 4","Random pseudofeature 5"),values = palette_OkabeIto_6)+
#   theme_bw(base_size = 12)+
#   xlab("Training epoch")+
#   ylab("L\u00b2 norm")+
#   guides(colour = guide_legend(title = "Feature set",override.aes = list(alpha = 1)))+
#   facet_wrap(~model_params,ncol = 3,labeller = facet_titles) +
#   theme(strip.background = element_rect(fill="white"))
# 
# ggsave("L2_norm_plot.png",height = 2440,width = 2880,unit="px")
# 
# ggplot(all_embedding_dist,aes(x = epoch,y = mantel_correlation,col=feature_set)) +
#   geom_point()+
#   geom_line()+
#   scale_colour_manual(labels = c("GC %","Dinucleotide composition","3mer composition","Amino acid composition",
#                                  "Dipeptide composition","Random pseudofeature"),values = palette_OkabeIto_6)+
#   theme_bw(base_size = 12)+
#   xlab("Training epoch")+
#   ylab("Mantel correlation")+
#   guides(colour = guide_legend(title = "Feature set",override.aes = list(alpha = 1)))+
#   facet_wrap(~model_params,ncol = 3,labeller = facet_titles) +
#   theme(strip.background = element_rect(fill="white"))
# 
# ggsave("mantel_correlation_plot.png",height = 2440,width = 2880,unit="px")

### now plot with accuracy/loss on a shared axis

facet_titles_feats <- labeller(feature_set = 
                                 c("dponly" = "Dipeptide composition",
                                   "rand" = "Random pseudofeature"))

all_training_logs_norand %<>%
  dplyr::rename("feature_set" = feat)

ggplot(all_embedding_dist_norand,aes(x = epoch,y = sum_square_residuals_norm,col=feature_set)) +
  # geom_point(size = 0.5)+
  geom_line(linewidth = 1)+
  geom_pointrange(data = rand_error,aes(x = epoch,y=sum_square_mean_norm,ymin = sum_square_mean_norm-sum_square_sd_norm,ymax = sum_square_mean_norm+sum_square_sd_norm),size = 0.02)+
  geom_line(data = rand_error,aes(x = epoch,y=sum_square_mean_norm),linewidth = 1)+
  # geom_point(data=all_training_logs,aes(x=epoch,y=test_acc,col=feature_set),shape = 4,size = 0.5)+
  geom_line(data=all_training_logs_norand,aes(x=epoch,y=all_acc_mean,col=feature_set),linetype = "dashed",linewidth = 1)+
  geom_pointrange(data = rand_training_logs,aes(x = epoch,y=all_acc_mean,ymin = all_acc_mean-all_acc_sd,ymax = all_acc_mean+all_acc_sd),size = 0.02,shape = 4)+
  geom_line(data = rand_training_logs,aes(x = epoch,y=all_acc_mean),linetype = "dashed",linewidth = 1)+
  scale_colour_manual(labels = c("GC %","Dinucleotide composition","3mer composition","Amino acid composition",
                                 "Dipeptide composition","Random pseudofeature"),values = palette_OkabeIto_2)+
  theme_bw(base_size = 12)+
  xlab("Training epoch")+
  ylab("Sum of squared residuals")+
  guides(colour = guide_legend(title = "Feature set",override.aes = list(alpha = 1))) +
  facet_wrap(~feature_set,ncol = 3,labeller = facet_titles_feats) +
  theme(strip.background = element_rect(fill="white"),panel.grid.minor = element_blank())

ggsave("./results_summary/sum_square_residuals_accuracy_plot_nopoint.png",height = 2000,width = 3500,unit="px")


# ggplot(all_embedding_dist,aes(x = epoch,y = mantel_correlation,col=feature_set)) +
#   geom_point()+
#   geom_line()+
#   geom_point(data=all_training_logs,aes(x=Epoch,y=test_acc,col=feat),shape = 4)+
#   geom_line(data=all_training_logs,aes(x=Epoch,y=test_acc,col=feat),linetype = "dashed")+
#   scale_colour_manual(labels = c("GC %","Dinucleotide composition","3mer composition","Amino acid composition",
#                                  "Dipeptide composition","Random pseudofeature 1","Random pseudofeature 2",
#                                  "Random pseudofeature 3","Random pseudofeature 4","Random pseudofeature 5"),values = palette_OkabeIto_6)+
#   theme_bw(base_size = 12)+
#   xlab("Training epoch")+
#   ylab("Mantel correlation")+
#   guides(colour = guide_legend(title = "Feature set")) +
#   facet_wrap(~model_params,ncol = 3,labeller = facet_titles) +
#   theme(strip.background = element_rect(fill="white"))
# 
# ggsave("mantel_accuracy_plot.png",height = 2440,width = 2880,unit="px")

# plot fewer points for clarity
# all_embedding_dist_norand_10 <- all_embedding_dist_norand %>%
#   filter(epoch %% 10 == T)

rand_error_10 <- rand_error %>%
  filter(epoch %% 10 == T)

# all_training_logs_norand_10 <- all_training_logs_norand %>%
#   filter(epoch %% 10 == T)

rand_training_logs_10 <- rand_training_logs %>%
  filter(epoch %% 10 == T)

ggplot(all_embedding_dist_norand,aes(x = epoch,y = sum_square_residuals_norm,col=feature_set)) +
  # geom_point(data = all_embedding_dist_10,aes(x = epoch,y = sum_square_residuals_norm,col=feature_set))+
  geom_line(linewidth = 0.8)+
  geom_pointrange(data = rand_error_10,aes(x = epoch,y=sum_square_mean_norm,ymin = sum_square_mean_norm-sum_square_sd_norm,ymax = sum_square_mean_norm+sum_square_sd_norm),size = 0.02)+
  geom_line(data = rand_error,aes(x = epoch,y=sum_square_mean_norm),linewidth = 0.8)+
  # geom_point(data=all_training_logs_10,aes(x=epoch,y=test_acc,col=feature_set),shape = 4)+
  geom_line(data=all_training_logs_norand,aes(x=epoch,y=all_acc_mean,col=feature_set),linetype = "dashed",linewidth = 0.8)+
  geom_pointrange(data = rand_training_logs_10,aes(x = epoch,y=all_acc_mean,ymin = all_acc_mean-all_acc_sd,ymax = all_acc_mean+all_acc_sd),size = 0.02,shape = 4)+
  geom_line(data = rand_training_logs,aes(x = epoch,y=all_acc_mean),linetype = "dashed",linewidth = 0.8)+
  scale_colour_manual(labels = c("GC %","Dinucleotide composition","3mer composition","Amino acid composition",
                                 "Dipeptide composition","Random pseudofeature"),values = palette_OkabeIto_2)+
  theme_bw(base_size = 12)+
  xlab("Training epoch")+
  ylab("Normalised sum of squared residuals / Accuracy (%)")+
  guides(colour = guide_legend(title = "Feature set")) +
  facet_wrap(~feature_set,labeller = facet_titles_feats) +
  theme(strip.background = element_rect(fill="white"),panel.grid.minor = element_blank())

ggsave("./results_summary/sum_square_residuals_accuracy_plot_randfewpoints.png",height = 2000,width = 3500,unit="px")


