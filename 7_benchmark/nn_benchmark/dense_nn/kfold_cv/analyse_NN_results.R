library(tidyverse)
library(readxl)
library(viridis)
library(magrittr)

`%nin%` <- Negate(`%in%`)

# read in subgraph results
agg_results <- read_xlsx(paste0("./results_summary/subgraph_agg_performance_all_params.xlsx"))
roc_results <- read_xlsx(paste0("./results_summary/subgraph_roc_performance_all_params.xlsx"))
noinfo_results <- read_xlsx(paste0("./results_summary/subgraph_noinfo_performance_all_params.xlsx"))

agg_results$feat <- factor(agg_results$feat,levels = c("dponly","rand1","rand2","rand3","rand4","rand5"))
roc_results$feat <- factor(roc_results$feat,levels = c("dponly","rand1","rand2","rand3","rand4","rand5"))
noinfo_results$feat <- factor(noinfo_results$feat,levels = c("dponly","rand1","rand2","rand3","rand4","rand5"))

process_rand <- function(x){
  rand_results <- x %>%
    filter(feat %in% c("rand1","rand2","rand3","rand4","rand5"))
  
  # take mean results across the five random pseudofeatures
  rand_results %<>%
    group_by(params,tax_level,Subgraph) %>%
    summarise(mean_accuracy_sd = sd(mean_accuracy),
              mean_accuracy = mean(mean_accuracy),
              mean_macro_sd = sd(mean_macro),
              mean_macro = mean(mean_macro))
  
  rand_results$feat <- "rand"
  
  feat_results <- x %>%
    filter(feat %nin% c("rand1","rand2","rand3","rand4","rand5"))
  
  feat_results$mean_accuracy_sd <- NA
  feat_results$mean_macro_sd <- NA
  
  feat_results <- rbind(feat_results,rand_results)
  return(feat_results)
}

agg_results <- process_rand(agg_results)
noinfo_results <- process_rand(noinfo_results)

# read in subgraph IDs
subgraph_id <- read_xlsx("../../../../5_phylogeny/3_generate_graph/subgraph_virus_classes.xlsx")

agg_results <- left_join(agg_results,subgraph_id)
roc_results <- left_join(roc_results,subgraph_id)
noinfo_results <- left_join(noinfo_results,subgraph_id)

subgraph_levels <- c("Whole Graph","Mammalian rubulaviruses","Mammalian morbilliviruses",
                     "Mammalian respiroviruses & fish viruses",
                     "Mammalian narmoviruses","Mammalian jeilongviruses",
                     "Mammalian henipaviruses","Avian orthoavulaviruses","Avian metaavulaviruses",
                     "Paraavulavirus wisconsinense","Reptile viruses")

agg_results$Members <- factor(agg_results$Members,levels = subgraph_levels)

roc_results$Members <- factor(roc_results$Members,levels = subgraph_levels)

noinfo_results$Members <- factor(noinfo_results$Members,levels = subgraph_levels)

palette <- RColorBrewer::brewer.pal(n = 6,"Set1")

facet_titles <- labeller(feat =
                           c("dponly" = "Dipeptide composition",
                             "rand" = "Random pseudofeature"))

agg_results %<>%
  separate_wider_delim(params,delim = "_",names = c("nhidden","nlayers"))

roc_results %<>%
  separate_wider_delim(params,delim = "_",names = c("nhidden","nlayers"))

noinfo_results %<>%
  separate_wider_delim(params,delim = "_",names = c("nhidden","nlayers"))

# factorise
hidden_levels <- c("32","64","128","256","512","1024")

agg_results$nhidden <- factor(agg_results$nhidden,levels = hidden_levels)

roc_results$nhidden <- factor(roc_results$nhidden,levels = hidden_levels)

noinfo_results$nhidden <- factor(noinfo_results$nhidden,levels = hidden_levels)

layers_levels <- c("2","3","4")

agg_results$nlayers <- factor(agg_results$nlayers,levels = layers_levels)

roc_results$nlayers <- factor(roc_results$nlayers,levels = layers_levels)

noinfo_results$nlayers <- factor(noinfo_results$nlayers,levels = layers_levels)

# split no info tables into order and class
noinfo_order <- noinfo_results %>% 
  dplyr::filter(tax_level == "Order") %>% 
  distinct()

noinfo_class <- noinfo_results %>% 
  dplyr::filter(tax_level == "Class") %>% 
  distinct()

ggplot(agg_results %>% dplyr::filter(tax_level == "Order"),aes(x = Members, y = mean_accuracy)) +
  geom_point(aes(col = nhidden,shape = nlayers),size = 2) +
  scale_color_manual(values = palette) +
  scale_shape_manual(values = c(21,22,24)) +
  geom_point(data = noinfo_order,aes(x = Members, y = mean_accuracy),col="grey3",shape = 4) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0,1,0.1))+
  theme(strip.background = element_rect(fill="white"),
        axis.text.x = element_text(angle = 315,hjust=0),
        panel.grid.minor = element_blank(),
        legend.position = "top") +
  guides(col=guide_legend(title="Dimension size"),shape = guide_legend(title="Hidden layers")) +
  facet_wrap(~feat,labeller = facet_titles)

ggsave("./results_summary/single_ablation_simple_feats_order_byparams.png",width = 2500,height=2500,units="px")

ggplot(agg_results %>% dplyr::filter(tax_level == "Order"),aes(x = feat, y = mean_accuracy)) +
  geom_point(aes(col = nhidden,shape = nlayers),size = 2) +
  scale_color_manual(values = palette) +
  scale_shape_manual(values = c(21,22,24)) +
  geom_hline(data = noinfo_order,aes(yintercept = mean_accuracy),col="grey3",linetype = "dashed") +
  theme_bw() +
  scale_y_continuous(breaks = seq(0,1,0.1))+
  theme(strip.background = element_rect(fill="white"),
        axis.text.x = element_text(angle = 315,hjust=0),
        panel.grid.minor = element_blank()) +
  guides(col=guide_legend(title="Dimension size"),shape = guide_legend(title="Hidden layers")) +
  facet_wrap(~Members,labeller = facet_titles)

ggsave("./results_summary/single_ablation_simple_feats_order_bysubgraph_byparams.png",width = 2500,height=2500,units="px")

acc_tbl <- agg_results %>% 
  filter(tax_level == "Order" & Members == "Whole Graph") %>% 
  group_by(nhidden,nlayers,feat) %>% 
  summarise(mean_acc = mean(mean_accuracy))

bt_tbl <- acc_tbl %>% 
  group_by(feat) %>% 
  filter(mean_acc == max(mean_acc))

bt_epochs <- tibble(feat = c(as.vector(bt_tbl$feat),"rand2","rand3","rand4","rand5"),epochs = NA)

bt_epochs$feat[bt_epochs$feat == "rand"] <- "rand1"

# read in best tune outputs to get training epochs
for(i in 1:nrow(bt_epochs)){
  
  # for the random variables average best # of training epochs across all random feature sets
  if(grepl("rand",bt_epochs$feat[i]) == T){
    rand_index <- bt_tbl$feat == "rand"
    bt_output <- read_csv(paste0("./feature_sets/results_",bt_epochs$feat[i],"/nn_output_",bt_tbl$nhidden[rand_index],"_",bt_tbl$nlayers[rand_index],".csv"))
  } else {bt_output <- read_csv(paste0("./feature_sets/results_",bt_tbl$feat[i],"/nn_output_",bt_tbl$nhidden[i],"_",bt_tbl$nlayers[i],".csv"))}
  
  bt_epochs$epochs[i] <- round(mean(bt_output$Training_epochs))
}

rand_mean_epochs <- bt_epochs %>% 
  filter(grepl("rand",feat)) %>% 
  summarise(epochs = round(mean(epochs)))

rand_mean_epochs$feat <- "rand"

bt_epochs %<>% 
  filter(!grepl("rand",feat))

bt_epochs <- rbind(bt_epochs,rand_mean_epochs)

bt_tbl <- left_join(bt_tbl,bt_epochs)

bt_tbl$nhidden <- as.numeric(as.character(bt_tbl$nhidden))

bt_tbl %<>%
  filter(case_when(feat == "dponly" ~ nhidden == min(nhidden) & epochs == min(epochs),
                   feat == "rand" ~ mean_acc == min(mean_acc)))

write_csv(bt_tbl,"./results_summary/best_tune_acc.csv")

max_tbl <- agg_results %>% 
  # filter(tax_level == "Order" & Members == "Whole Graph") %>% 
  group_by(feat,tax_level,Members) %>% 
  summarise(max_acc = max(mean_accuracy))

write_csv(max_tbl,"./results_summary/max_accuracy_full.csv")

# filter for the best performing tunes on the whole graph
whole_graph_best_tune <- agg_results %>%
  group_by(feat,tax_level,Members) %>% 
  filter(mean_accuracy == max(mean_accuracy)) %>% 
  filter(tax_level == "Order" & Members == "Whole Graph") %>% 
  ungroup()

dp_bt <- whole_graph_best_tune %>% filter(feat == "dponly")
rand_bt <- whole_graph_best_tune %>% filter(feat == "rand")

plot_tbl <- agg_results %>% 
  filter(feat == "dponly" & tax_level == "Order" & nhidden == dp_bt$nhidden & nlayers == dp_bt$nlayers |
           feat == "rand" & tax_level == "Order" & nhidden == rand_bt$nhidden & nlayers == rand_bt$nlayers)

# plot only the best tune -------------------------------------------------
palette_OkabeIto_6 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                        "#CC79A7", "#999999")

plot_tbl %<>%
  mutate(feat = case_when(feat == "dponly" ~ "Dipeptide composition",
                          feat == "rand" ~ "Random pseudofeature"))

plot_tbl$feat <- factor(plot_tbl$feat,levels = c("Dipeptide composition","Random pseudofeature"))

benchmark_tbl <- tibble(accuracy = 0.9255814,Members = "Whole Graph")

plot_tbl_feats <- plot_tbl %>%
  filter(feat != "Random pseudofeature")

plot_tbl_rand <- plot_tbl %>%
  filter(feat == "Random pseudofeature")

ggplot(plot_tbl_feats,aes(x = Members, y = mean_accuracy,col=feat)) +
  geom_point(size = 3.5) +
  geom_pointrange(data = plot_tbl_rand, aes(x = Members, y = mean_accuracy,col=feat,ymin = mean_accuracy - mean_accuracy_sd,ymax = mean_accuracy + mean_accuracy_sd),size = 0.8)+
  scale_color_manual(values = palette_OkabeIto_6) +
  # geom_point(data = noinfo_order,aes(x = Members, y = mean_accuracy),col="grey3",shape = 4,size = 3) +
  geom_point(data = benchmark_tbl,aes(x = Members, y = accuracy),col="magenta2",shape = 8,size = 3) +
  theme_bw() +
  ylab("Accuracy")+
  scale_y_continuous(breaks = seq(0.5,1,0.1),limits = c(0.5,1))+
  theme(strip.background = element_rect(fill="white"),
        axis.text.x = element_text(angle = 315,hjust=0),
        panel.grid.minor = element_blank()) +
  guides(col=guide_legend(title="Feature set"))

ggsave("./results_summary/best_tune_subgraph_comparison_no_noinfo.png",width = 2500,height=2500,units="px")
