library(tidyverse)
library(readxl)
library(viridis)
library(magrittr)

`%nin%` <- Negate(`%in%`)

# read in subgraph results
agg_results <- read_xlsx(paste0("./results_summary/subgraph_agg_performance_all_params.xlsx"))
roc_results <- read_xlsx(paste0("./results_summary/subgraph_roc_performance_all_params.xlsx"))
noinfo_results <- read_xlsx(paste0("./results_summary/subgraph_noinfo_performance_all_params.xlsx"))

agg_results$feat <- factor(agg_results$feat,levels = c("gconly","dintonly","kmer3only","aaconly","dponly","rand1","rand2","rand3","rand4","rand5"))
roc_results$feat <- factor(roc_results$feat,levels = c("gconly","dintonly","kmer3only","aaconly","dponly","rand1","rand2","rand3","rand4","rand5"))
noinfo_results$feat <- factor(noinfo_results$feat,levels = c("gconly","dintonly","kmer3only","aaconly","dponly","rand1","rand2","rand3","rand4","rand5"))

process_rand <- function(x){
  rand_results <- x %>%
    filter(feat %in% c("rand1","rand2","rand3","rand4","rand5"))
  
  # take mean results across the five random pseudofeatures
  rand_results %<>%
    group_by(dropout,tax_level,Subgraph) %>%
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

subgraph_levels <- c("Whole Graph","Mammalian rubulaviruses","Mammalian henipaviruses",
                     "Mammalian respiroviruses & fish viruses", "Mammalian morbilliviruses",
                     "Mammalian jeilongviruses","Mammalian narmoviruses",
                     "Avian orthoavulaviruses","Avian metaavulaviruses",
                     "Reptile viruses","Paraavulavirus wisconsinense")

agg_results$Members <- factor(agg_results$Members,levels = subgraph_levels)

roc_results$Members <- factor(roc_results$Members,levels = subgraph_levels)

noinfo_results$Members <- factor(noinfo_results$Members,levels = subgraph_levels)

palette <- RColorBrewer::brewer.pal(n = 4,"Set1")

facet_titles <- labeller(feat =
                           c("gconly" = "GC %",
                             "aaconly" = "Amino acid composition",
                             "dintonly" = "Dinucleotide composition",
                             "dponly" = "Dipeptide composition",
                             "kmer3only" = "3mers",
                             "rand" = "Random pseudofeature"))

# factorise
dropout_levels <- c("0.2","0.3","0.4","0.5")

agg_results$dropout <- factor(agg_results$dropout,levels = dropout_levels)

roc_results$dropout <- factor(roc_results$dropout,levels = dropout_levels)

noinfo_results$dropout <- factor(noinfo_results$dropout,levels = dropout_levels)

# split no info tables into order and class
noinfo_order <- noinfo_results %>% 
  dplyr::filter(tax_level == "Order") %>% 
  distinct()

noinfo_class <- noinfo_results %>% 
  dplyr::filter(tax_level == "Class") %>% 
  distinct()

ggplot(agg_results %>% dplyr::filter(tax_level == "Order"),aes(x = Members, y = mean_accuracy)) +
  geom_point(aes(col = dropout),size = 2) +
  scale_color_manual(values = palette) +
  geom_point(data = noinfo_order,aes(x = Members, y = mean_accuracy),col="grey3",shape = 4) +
  theme_bw() +
  ylab("Accuracy")+
  xlab(NULL)+
  scale_y_continuous(breaks = seq(0,1,0.1))+
  theme(strip.background = element_rect(fill="white"),
        axis.text.x = element_text(angle = 315,hjust=0),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        plot.margin = margin(0,2.5,0,0, "cm")) +
  guides(col=guide_legend(title="Dropout rate")) +
  facet_wrap(~feat,labeller = facet_titles)

ggsave("./results_summary/single_ablation_simple_feats_order_bydropout.png",width = 2500,height=2500,units="px")

ggplot(agg_results %>% dplyr::filter(tax_level == "Order"),aes(x = Members, y = mean_accuracy)) +
  geom_point(aes(col = dropout),size = 2) +
  scale_color_manual(values = palette) +
  theme_bw() +
  ylab("Accuracy")+
  xlab(NULL)+
  scale_y_continuous(breaks = seq(0,1,0.1))+
  theme(strip.background = element_rect(fill="white"),
        axis.text.x = element_text(angle = 315,hjust=0),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        plot.margin = margin(0,2.5,0,0, "cm")) +
  guides(col=guide_legend(title="Dropout rate")) +
  facet_wrap(~feat,labeller = facet_titles)

ggsave("./results_summary/single_ablation_simple_feats_order_bydropout_nonoinfo.png",width = 2500,height=2500,units="px")

ggplot(agg_results %>% dplyr::filter(tax_level == "Order"),aes(x = feat, y = mean_accuracy)) +
  geom_point(aes(col = dropout),size = 2) +
  scale_color_manual(values = palette) +
  geom_hline(data = noinfo_order,aes(yintercept = mean_accuracy),col="grey3",linetype = "dashed") +
  theme_bw() +
  ylab("Accuracy")+
  scale_y_continuous(breaks = seq(0,1,0.1))+
  theme(strip.background = element_rect(fill="white"),
        axis.text.x = element_text(angle = 315,hjust=0),
        panel.grid.minor = element_blank()) +
  guides(col=guide_legend(title="Dropout rate")) +
  facet_wrap(~Members,labeller = facet_titles) 


ggsave("./results_summary/single_ablation_simple_feats_order_bysubgraph_bydropout.png",width = 2500,height=2500,units="px")

acc_tbl <- agg_results %>% 
  filter(tax_level == "Order" & Members == "Whole Graph") %>% 
  group_by(dropout,feat) %>% 
  summarise(mean_acc = mean(mean_accuracy))

write_csv(acc_tbl %>% 
            group_by(feat) %>% 
            filter(mean_acc == max(mean_acc)),"./results_summary/best_tune_acc.csv")

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

aac_bt <- whole_graph_best_tune %>% filter(feat == "aaconly")
dint_bt <- whole_graph_best_tune %>% filter(feat == "dintonly")
dp_bt <- whole_graph_best_tune %>% filter(feat == "dponly")
gc_bt <- whole_graph_best_tune %>% filter(feat == "gconly")
kmer3_bt <- whole_graph_best_tune %>% filter(feat == "kmer3only")
rand_bt <- whole_graph_best_tune %>% filter(feat == "rand")

plot_tbl <- agg_results %>% 
  filter(feat == "aaconly" & tax_level == "Order" & dropout == aac_bt$dropout |
          feat == "dintonly" & tax_level == "Order" & dropout == dint_bt$dropout |
           feat == "dponly" & tax_level == "Order" & dropout == dp_bt$dropout |
           feat == "gconly" & tax_level == "Order" & dropout == gc_bt$dropout |
           feat == "kmer3only" & tax_level == "Order" & dropout == kmer3_bt$dropout |
           feat == "rand" & tax_level == "Order" & dropout == rand_bt$dropout)

# plot only the best tune -------------------------------------------------
palette_OkabeIto_6 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                        "#CC79A7", "#999999")

plot_tbl %<>%
  mutate(feat = case_when(feat == "dintonly" ~ "Dinucleotide composition",
                          feat == "aaconly" ~ "Amino acid composition",
                          feat == "dponly" ~ "Dipeptide composition",
                          feat == "kmer3only" ~ "3mer composition",
                          feat == "gconly" ~ "GC %",
                          feat == "rand" ~ "Random pseudofeature"))

plot_tbl$feat <- factor(plot_tbl$feat,levels = c("GC %","Dinucleotide composition","3mer composition",
                                                       "Amino acid composition","Dipeptide composition","Random pseudofeature"))

benchmark_tbl <- tibble(accuracy = 0.9255814,Members = "Whole Graph")

plot_tbl_feats <- plot_tbl %>%
  filter(feat != "Random pseudofeature")

plot_tbl_rand <- plot_tbl %>%
  filter(feat == "Random pseudofeature")

# in some cases multiple tunes have the same results - filter these out
view(plot_tbl_feats %>% 
       dplyr::count(feat,Subgraph))

plot_tbl_feats %<>% 
  filter(case_when(feat == "Dipeptide composition" ~ dropout == 0.2,
                   T ~ dropout %in% unique(plot_tbl_feats$dropout))) # this accepts all observations that do not match previous condition

plot_tbl_feats %<>% 
  filter(Members != "Reptile viruses" & Members != "Paraavulavirus wisconsinense")

plot_tbl_rand %<>%
  filter(Members != "Reptile viruses" & Members != "Paraavulavirus wisconsinense")

ggplot(plot_tbl_feats,aes(x = Members, y = mean_accuracy,col=feat)) +
  geom_point(size = 4,position = position_dodge(width = 0.4)) +
  geom_pointrange(data = plot_tbl_rand, aes(x = Members, y = mean_accuracy,col=feat,ymin = mean_accuracy - mean_accuracy_sd,ymax = mean_accuracy + mean_accuracy_sd),size = 0.75)+
  scale_color_manual(values = palette_OkabeIto_6) +
  # geom_point(data = noinfo_order,aes(x = Members, y = mean_accuracy),col="grey3",shape = 4,size = 3) +
  geom_point(data = benchmark_tbl,aes(x = Members, y = accuracy),col="black",shape = 8,size = 4) +
  theme_bw(base_size = 15) +
  ylab("Accuracy")+
  xlab(NULL)+
  scale_y_continuous(breaks = seq(0.5,1,0.1),limits = c(0.5,1))+
  theme(strip.background = element_rect(fill="white"),
        axis.text.x = element_text(angle = 315,hjust=0),
        legend.position = "top",
        panel.grid.minor = element_blank(),
        plot.margin = margin(0,2.5,0,0.2, "cm")) +
  guides(col=guide_legend(title="Feature set"))

ggsave("./results_summary/best_tune_subgraph_comparison_no_noinfo_short.pdf",width = 5000,height=5500,units="px",dpi=600)
