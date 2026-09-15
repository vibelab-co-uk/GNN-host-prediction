library(tidyverse)
library(readxl)
library(viridis)
library(magrittr)

`%nin%` <- Negate(`%in%`)

# read in subgraph results
agg_results <- read_xlsx(paste0("./results_summary/subgraph_agg_performance_all_params_trailer.xlsx"))
roc_results <- read_xlsx(paste0("./results_summary/subgraph_roc_performance_all_params_trailer.xlsx"))
noinfo_results <- read_xlsx(paste0("./results_summary/subgraph_noinfo_performance_all_params_trailer.xlsx"))

agg_results$feat <- factor(agg_results$feat,levels = c("gconly","dintonly","kmer3only","kmer4only","kmer5only","kmer6only"))
roc_results$feat <- factor(roc_results$feat,levels = c("gconly","dintonly","kmer3only","kmer4only","kmer5only","kmer6only"))
noinfo_results$feat <- factor(noinfo_results$feat,levels = c("gconly","dintonly","kmer3only","kmer4only","kmer5only","kmer6only"))

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

# facet_titles <- labeller(feat =
#                            c("gconly" = "GC %",
#                              "aaconly" = "Amino acid composition",
#                              "dintonly" = "Dinucleotide composition",
#                              "dponly" = "Dipeptide composition",
#                              "kmer3only" = "3mers",
#                              "rand" = "Random pseudofeature"),
#                          params = 
#                            c("128_2" = "Layers: 2 Hidden dims: 128",
#                              "128_3" = "Layers: 3 Hidden dims: 128",
#                              "128_4" = "Layers: 4 Hidden dims: 128",
#                              "256_2" = "Layers: 2 Hidden dims: 256",
#                              "256_3" = "Layers: 3 Hidden dims: 256",
#                              "256_4" = "Layers: 4 Hidden dims: 256",
#                              "512_2" = "Layers: 2 Hidden dims: 512",
#                              "512_3" = "Layers: 3 Hidden dims: 512",
#                              "512_4" = "Layers: 4 Hidden dims: 512"))

facet_titles <- labeller(feat =
                           c("gconly" = "GC %",
                             "dintonly" = "Dinucleotide composition",
                             "kmer3only" = "3mers",
                             "kmer4only" = "4mers",
                             "kmer5only" = "5mers",
                             "kmer6only" = "6mers"))

agg_results %<>%
  separate_wider_delim(params,delim = "_",names = c("nhidden","nlayers"))

roc_results %<>%
  separate_wider_delim(params,delim = "_",names = c("nhidden","nlayers"))

noinfo_results %<>%
  separate_wider_delim(params,delim = "_",names = c("nhidden","nlayers"))

# factorise
hidden_levels <- c("128","256","512","1024")

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

ggsave("./results_summary/single_ablation_simple_feats_order_byparams_trailerUTR.png",width = 2500,height=2500,units="px")

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

ggsave("./results_summary/single_ablation_simple_feats_order_bysubgraph_byparams_trailerUTR.png",width = 2500,height=2500,units="px")

acc_tbl <- agg_results %>% 
  filter(tax_level == "Order" & Members == "Whole Graph") %>% 
  group_by(nhidden,nlayers,feat) %>% 
  summarise(mean_acc = mean(mean_accuracy))

write_csv(acc_tbl %>% 
            group_by(feat) %>% 
            filter(mean_acc == max(mean_acc)),"./results_summary/best_tune_acc_trailerUTR.csv")

max_tbl <- agg_results %>% 
  # filter(tax_level == "Order" & Members == "Whole Graph") %>% 
  group_by(feat,tax_level,Members) %>% 
  summarise(max_acc = max(mean_accuracy))

write_csv(max_tbl,"./results_summary/max_accuracy_full_trailerUTR.csv")

# filter for the best performing tunes on the whole graph
whole_graph_best_tune <- agg_results %>%
  group_by(feat,tax_level,Members) %>% 
  filter(mean_accuracy == max(mean_accuracy)) %>% 
  filter(tax_level == "Order" & Members == "Whole Graph") %>% 
  ungroup()

dint_bt <- whole_graph_best_tune %>% filter(feat == "dintonly")
gc_bt <- whole_graph_best_tune %>% filter(feat == "gconly")
kmer3_bt <- whole_graph_best_tune %>% filter(feat == "kmer3only")
kmer4_bt <- whole_graph_best_tune %>% filter(feat == "kmer4only")
kmer5_bt <- whole_graph_best_tune %>% filter(feat == "kmer5only")
kmer6_bt <- whole_graph_best_tune %>% filter(feat == "kmer6only")

# # where there are multiple tunes with identical performance, choose the lower dimensional model
dint_bt %<>%
  filter(nhidden == "512")

kmer4_bt %<>%
  filter(nhidden == "512")

kmer6_bt %<>%
  filter(nhidden == "128")


plot_tbl <- agg_results %>% 
  filter(feat == "dintonly" & tax_level == "Order" & nhidden == dint_bt$nhidden & nlayers == dint_bt$nlayers |
           feat == "gconly" & tax_level == "Order" & nhidden == gc_bt$nhidden & nlayers == gc_bt$nlayers|
           feat == "kmer3only" & tax_level == "Order" & nhidden == kmer3_bt$nhidden & nlayers == kmer3_bt$nlayers |
           feat == "kmer4only" & tax_level == "Order" & nhidden == kmer4_bt$nhidden & nlayers == kmer4_bt$nlayers |
           feat == "kmer5only" & tax_level == "Order" & nhidden == kmer5_bt$nhidden & nlayers == kmer5_bt$nlayers |
           feat == "kmer6only" & tax_level == "Order" & nhidden == kmer6_bt$nhidden & nlayers == kmer6_bt$nlayers)

# plot only the best tune -------------------------------------------------
palette_OkabeIto_6 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                        "#CC79A7", "#0072B2")

plot_tbl %<>%
  mutate(feat = case_when(feat == "dintonly" ~ "Dinucleotide composition",
                          feat == "gconly" ~ "GC %",
                          feat == "kmer3only" ~ "3mer composition",
                          feat == "kmer4only" ~ "4mer composition",
                          feat == "kmer5only" ~ "5mer composition",
                          feat == "kmer6only" ~ "6mer composition"))

plot_tbl$feat <- factor(plot_tbl$feat,levels = c("GC %","Dinucleotide composition","3mer composition",
                                                 "4mer composition","5mer composition","6mer composition"))

benchmark_tbl <- tibble(accuracy = 0.9255814,Members = "Whole Graph")

ggplot(plot_tbl,aes(x = Members, y = mean_accuracy,col=feat)) +
  geom_point(size = 3.5) +
  scale_color_manual(values = palette_OkabeIto_6) +
  # geom_point(data = noinfo_order,aes(x = Members, y = mean_accuracy),col="grey3",shape = 4,size = 3) +
  geom_point(data = benchmark_tbl,aes(x = Members, y = accuracy),col="magenta2",shape = 8,size = 3) +
  theme_bw() +
  ylab("Accuracy")+
  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(0,1))+
  theme(strip.background = element_rect(fill="white"),
        axis.text.x = element_text(angle = 315,hjust=0),
        panel.grid.minor = element_blank()) +
  guides(col=guide_legend(title="Feature set"))

ggsave("./results_summary/best_tune_subgraph_comparison_no_noinfo_trailerUTR.png",width = 2500,height=2500,units="px")

