library(tidyverse)
library(readxl)
library(viridis)
library(magrittr)
library(fishualize)

`%nin%` <- Negate(`%in%`)

### plot results of reintro blocked cross validation for dense NN and GNN with BLAST benchmark

# read in subgraph results
agg_results_nn <- read_xlsx(paste0("./results_summary/subgraph_agg_performance_best_tune.xlsx"))
agg_results_nn_0 <- read_xlsx(paste0("../blocked_cv/results_summary/subgraph_agg_performance_best_tune_testonly.xlsx"))

agg_results_gnn <- read_xlsx(paste0("../../../../6_GNN/3_graph_modelling/blocked_cv_resample/simple_feats/results_summary/subgraph_agg_performance_best_tune.xlsx"))
agg_results_gnn_0 <- read_xlsx(paste0("../../../../6_GNN/3_graph_modelling/blocked_cv/simple_feats/results_summary/subgraph_agg_performance_best_tune_testonly.xlsx"))

# filter for only dp + rand features
agg_results_gnn %<>%
  filter(feat %in% c("dponly","rand1","rand2","rand3","rand4","rand5"))

agg_results_gnn_0 %<>%
  filter(feat %in% c("dponly","rand1","rand2","rand3","rand4","rand5"))

# filter for the dropout value used to train resample models
agg_results_nn_0 %<>%
  filter(dropout == 0.2) %>% 
  select(-dropout)

agg_results_nn$model <- "FNN"
agg_results_nn_0$model <- "FNN"

agg_results_gnn$model <- "GNN"
agg_results_gnn_0$model <- "GNN"

agg_results <- rbind(agg_results_nn,agg_results_gnn)
agg_results_0 <- rbind(agg_results_nn_0,agg_results_gnn_0)

# read in blast benchmark results 
blast_results <- read_csv("../../../../7_benchmark/blastn_benchmark/reintro_blocked_cv/blast_reintro_accuracy_metrics.csv")
blast_results_0 <- read_csv("../../../../7_benchmark/blastn_benchmark/blocked_cv/blast_blocked_accuracy_metrics.csv")

blast_results_0 %<>%
  filter(Subgraph != 5 & Subgraph != 7)

# add in redundant columns for binding with other results
blast_results_0$N_samples <- 0
blast_results_0$sd_accuracy <- NA

blast_results_0 %<>%
  dplyr::rename("mean_accuracy" = "Accuracy")

blast_results <- rbind(blast_results,blast_results_0)

# process random repeats
results_rand <- agg_results %>% 
  filter(feat %in% c("rand1","rand2","rand3","rand4","rand5"))

agg_results <- agg_results %>% 
  filter(feat %nin% c("rand1","rand2","rand3","rand4","rand5"))

results_rand %<>% 
  group_by(Block,tax_level,N_samples,model) %>% 
  summarise(mean_macro = mean(mean_macro,na.rm=T),mean_accuracy = mean(mean_accuracy),sd_accuracy = mean(sd_accuracy))

results_rand$feat <- "rand"

agg_results <- rbind(agg_results,results_rand)

# process random repeats for no reintros
results_0_rand <- agg_results_0 %>% 
  filter(feat %in% c("rand1","rand2","rand3","rand4","rand5"))

agg_results_0 <- agg_results_0 %>% 
  filter(feat %nin% c("rand1","rand2","rand3","rand4","rand5"))

results_0_rand %<>% 
  group_by(Subgraph,tax_level,model) %>% 
  summarise(mean_macro = mean(mean_macro,na.rm=T),mean_accuracy = mean(mean_accuracy))
  
results_0_rand$feat <- "rand"

agg_results_0 <- rbind(agg_results_0,results_0_rand)

agg_results$feat <- factor(agg_results$feat,levels = c("gconly","dintonly","kmer3only","aaconly","dponly","rand"))
agg_results_0$feat <- factor(agg_results_0$feat,levels = c("gconly","dintonly","kmer3only","aaconly","dponly","rand"))

agg_results %<>% 
  dplyr::rename("Subgraph" = "Block")



subgraph_id <- read_xlsx("../../../../5_phylogeny/3_generate_graph/subgraph_virus_classes.xlsx")

agg_results <- left_join(agg_results,subgraph_id)
agg_results_0 <- left_join(agg_results_0,subgraph_id)
blast_results <- left_join(blast_results,subgraph_id)

agg_results_0$sd_accuracy <- NA
agg_results_0$N_samples <- 0

agg_results <- rbind(agg_results_0,agg_results)

subgraph_levels <- c("Mammalian rubulaviruses","Mammalian henipaviruses",
                     "Mammalian respiroviruses & fish viruses","Mammalian morbilliviruses",
                     "Mammalian jeilongviruses","Mammalian narmoviruses",
                     "Avian orthoavulaviruses","Avian metaavulaviruses")

agg_results$Members <- factor(agg_results$Members,levels = subgraph_levels)


palette_OkabeIto <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                      "#0072B2", "#D55E00", "#CC79A7", "#999999")


palette_OkabeIto_6 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                        "#CC79A7", "#999999")

agg_results %<>%
  mutate(feat = case_when(feat == "dponly" ~ "Dipeptide composition",
                            feat == "rand" ~ "Random pseudofeature"))

agg_results %<>% 
  filter(tax_level == "Order")

# plot with the blast benchmark
palette_OkabeIto_3 <- c("#E69F00", "#56B4E9","#999999")

palette_fish <- fish(n = 5, option = "Scarus_tricolor")
palette <- c(palette_fish[2],palette_fish[4],"#999999")

blast_results$feat <- "BLAST benchmark"

agg_results %<>%
  select(mean_accuracy,sd_accuracy,feat,Subgraph,N_samples,Members,model)

blast_results$model <- "BLAST benchmark"

agg_results <- rbind(agg_results,blast_results)

agg_results$feat <- factor(agg_results$feat,levels = c("Dipeptide composition","Random pseudofeature","BLAST benchmark"))

agg_results$model <- factor(agg_results$model,levels = c("GNN","FNN","BLAST benchmark"))

myshape <- c(1,2,4)

ggplot(agg_results,aes(x = N_samples, y = mean_accuracy,col = model,shape = feat)) +
  # geom_point(size = 2) +
  geom_line()+
  geom_pointrange(data = agg_results,
                  aes(x = N_samples,y=mean_accuracy,ymin = mean_accuracy-sd_accuracy,ymax = mean_accuracy+sd_accuracy),
                  size = 0.5)+
  scale_color_manual(values = palette) +
  scale_shape_manual(values = myshape)+
  theme_bw() +
  xlab("In-block train samples")+
  ylab("Mean accuracy")+
  scale_y_continuous(limits = c(0,1))+
  theme(strip.background = element_rect(fill="white"),
        panel.grid.minor = element_blank())+
  # legend.position = "top") +
  guides(col=guide_legend(title="Model type"),shape = guide_legend(title = "Feature")) +
  facet_wrap(~Members,nrow=2)

ggsave("./results_summary/reintro_accuracy_gain_wBLAST.png",width = 3500,height=1500,units="px")

myshape <- c(1,4)

ggplot(agg_results %>% filter(feat != "Random pseudofeature"),aes(x = N_samples, y = mean_accuracy,col = model,shape = feat)) +
  # geom_point(size = 2) +
  geom_line()+
  geom_pointrange(data = agg_results %>% filter(feat != "Random pseudofeature"),
                  aes(x = N_samples,y=mean_accuracy,ymin = mean_accuracy-sd_accuracy,ymax = mean_accuracy+sd_accuracy),
                  size = 0.5,linetype = "dotdash")+
  scale_color_manual(values = palette) +
  scale_shape_manual(values = myshape)+
  theme_bw() +
  xlab("In-block train samples")+
  ylab("Mean accuracy")+
  scale_y_continuous(limits = c(0,1))+
  theme(strip.background = element_rect(fill="white"),
        panel.grid.minor = element_blank())+
  # legend.position = "top") +
  guides(col=guide_legend(title="Model type"),shape = guide_legend(title = "Feature")) +
  facet_wrap(~Members,nrow=2)

ggsave("./results_summary/reintro_accuracy_gain_wBLAST_norand.png",width = 3500,height=1500,units="px")
