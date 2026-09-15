library(tidyverse)
library(readxl)
library(viridis)
library(magrittr)

`%nin%` <- Negate(`%in%`)

# read in subgraph results
agg_results <- read_xlsx(paste0("./results_summary/subgraph_agg_performance_best_tune_testonly.xlsx"))
roc_results <- read_xlsx(paste0("./results_summary/subgraph_roc_performance_best_tune_testonly.xlsx"))
noinfo_results <- read_xlsx(paste0("./results_summary/subgraph_noinfo_performance_best_tune_testonly.xlsx"))

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

subgraph_levels <- c("Whole Graph","Mammalian rubulaviruses","Mammalian morbilliviruses",
                     "Mammalian respiroviruses & fish viruses",
                     "Mammalian narmoviruses","Mammalian jeilongviruses",
                     "Mammalian henipaviruses","Avian orthoavulaviruses","Avian metaavulaviruses",
                     "Paraavulavirus wisconsinense","Reptile viruses")

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
  scale_y_continuous(breaks = seq(0,1,0.1))+
  theme(strip.background = element_rect(fill="white"),
        axis.text.x = element_text(angle = 315,hjust=0),
        panel.grid.minor = element_blank(),
        legend.position = "top") +
  guides(col=guide_legend(title="Dropout rate")) +
  facet_wrap(~feat,labeller = facet_titles)

ggsave("./results_summary/single_ablation_simple_feats_order_bydropout.png",width = 2500,height=2500,units="px")

ggplot(agg_results %>% dplyr::filter(tax_level == "Order"),aes(x = feat, y = mean_accuracy)) +
  geom_point(aes(col = dropout),size = 2) +
  scale_color_manual(values = palette) +
  geom_hline(data = noinfo_order,aes(yintercept = mean_accuracy),col="grey3",linetype = "dashed") +
  theme_bw() +
  scale_y_continuous(breaks = seq(0,1,0.1))+
  theme(strip.background = element_rect(fill="white"),
        axis.text.x = element_text(angle = 315,hjust=0),
        panel.grid.minor = element_blank()) +
  guides(col=guide_legend(title="Dropout rate")) +
  facet_wrap(~Members,labeller = facet_titles)

ggsave("./results_summary/single_ablation_simple_feats_order_bysubgraph_bydropout.png",width = 2500,height=2500,units="px")
