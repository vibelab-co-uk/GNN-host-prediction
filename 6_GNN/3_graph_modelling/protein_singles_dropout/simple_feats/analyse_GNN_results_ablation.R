library(tidyverse)
library(readxl)
library(viridis)
library(magrittr)
library(colorblindcheck)

`%nin%` <- Negate(`%in%`)

# read in subgraph results
agg_results <- read_xlsx("./results_summary/subgraph_agg_performance.xlsx")
# roc_results <- read_xlsx("./results_summary/subgraph_roc_performance.xlsx")
noinfo_results <- read_xlsx("./results_summary/subgraph_noinfo_performance.xlsx")

agg_results$feat <- factor(agg_results$feat,levels = c("gconly","dintonly","kmer3only","aaconly","dponly","rand1","rand2","rand3","rand4","rand5"))
# roc_results$feat <- factor(roc_results$feat,levels = c("gconly","dintonly","kmer3only","aaconly","dponly","rand1","rand2","rand3","rand4","rand5"))
noinfo_results$feat <- factor(noinfo_results$feat,levels = c("gconly","dintonly","kmer3only","aaconly","dponly","rand1","rand2","rand3","rand4","rand5"))

# take average of random sets
rand_results <- agg_results %>%
  filter(feat %in% c("rand1","rand2","rand3","rand4","rand5"))

rand_results %<>% 
  group_by(Ablation,tax_level,Subgraph) %>% 
  summarise(mean_accuracy = mean(mean_accuracy),mean_macro = mean(mean_macro))

rand_results$feat <- "rand"

agg_results %<>%
  filter(feat %nin% c("rand1","rand2","rand3","rand4","rand5"))

agg_results <- rbind(agg_results,rand_results)

# same for noinfo
rand_noinfo <- noinfo_results %>%
  filter(feat %in% c("rand1","rand2","rand3","rand4","rand5"))

rand_noinfo %<>% 
  group_by(Ablation,tax_level,Subgraph) %>% 
  summarise(mean_accuracy = mean(mean_accuracy),mean_macro = mean(mean_macro))

rand_noinfo$feat <- "rand"

noinfo_results %<>%
  filter(feat %nin% c("rand1","rand2","rand3","rand4","rand5"))

noinfo_results <- rbind(noinfo_results,rand_noinfo)

# read in full feature set results
agg_results_full <- read_csv("../../kfold_cv_early_stopping/simple_feats/results_summary/max_accuracy_full.csv")

agg_results_full$feat <- factor(agg_results_full$feat,levels = c("gconly","dintonly","kmer3only","aaconly","dponly","rand"))

subgraph_id <- read_xlsx("/Users/jamieherzig/Documents/clean_run/5_phylogeny/3_generate_graph/subgraph_virus_classes.xlsx")

agg_results <- left_join(agg_results,subgraph_id)
# roc_results <- left_join(roc_results,subgraph_id)
noinfo_results <- left_join(noinfo_results,subgraph_id)

agg_results_full <- left_join(agg_results_full,subgraph_id)

subgraph_levels <- c("Whole Graph","Mammalian rubulaviruses","Mammalian henipaviruses",
                     "Mammalian respiroviruses & fish viruses", "Mammalian morbilliviruses",
                     "Mammalian jeilongviruses","Mammalian narmoviruses",
                     "Avian orthoavulaviruses","Avian metaavulaviruses",
                     "Reptile viruses","Paraavulavirus wisconsinense")

agg_results$Members <- factor(agg_results$Members,levels = subgraph_levels)

# roc_results$Members <- factor(roc_results$Members,levels = subgraph_levels)

noinfo_results$Members <- factor(noinfo_results$Members,levels = subgraph_levels)

agg_results_full$Members <- factor(agg_results_full$Members,levels = subgraph_levels)

# roc_results_full$Members <- factor(roc_results_full$Members,levels = subgraph_levels)


facet_titles <- labeller(feat = 
                           c("gconly" = "GC %",
                             "aaconly" = "Amino acid composition",
                             "dintonly" = "Dinucleotide composition",
                             "dponly" = "Dipeptide composition",
                             "kmer3only" = "3mers",
                             "rand" = "Random pseudofeature"))

noinfo_order <- noinfo_results %>% 
  dplyr::filter(tax_level == "Order") %>% 
  dplyr::select(-Ablation) %>% 
  distinct()

noinfo_class <- noinfo_results %>% 
  dplyr::filter(tax_level == "Class") %>% 
  dplyr::select(-Ablation) %>% 
  distinct()

palette_OkabeIto <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                      "#0072B2", "#D55E00", "#CC79A7", "#999999")


palette_OkabeIto_6 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                        "#CC79A7", "#D55E00")

palette_OkabeIto_5 <- c("#56B4E9", "#009E73", "#F0E442", 
                        "#CC79A7", "#D55E00")

agg_results$Ablation <- factor(agg_results$Ablation,levels = c("Genomic","Attachment","Fusion","Nucleocapsid","Matrix","Polymerase"))

ggplot(agg_results %>% dplyr::filter(tax_level == "Order"),aes(x = Members, y = mean_accuracy)) +
  # geom_point(aes(col = Ablation),position = "jitter",size = 2) +
  geom_point(aes(col = Ablation),size = 2,shape = 1) +
  geom_point(data = noinfo_order,aes(x = Members, y = mean_accuracy),col="grey3",shape = 4) +
  geom_point(data = agg_results_full %>% dplyr::filter(tax_level == "Order"),aes(x = Members,y = max_acc),
             col = "magenta2",shape = 8,size = 1.5)+
  scale_color_manual(values = palette_OkabeIto_6) + 
  theme_bw() +
  scale_y_continuous(breaks = seq(0,1,0.1))+
  # guides(col = "none") +
  labs(col = "Feature source")+
  theme(strip.background = element_rect(fill="white"),
        axis.text.x = element_text(angle = 315,hjust=0),
        panel.grid.minor = element_blank(),
        legend.position = "top") +
  facet_wrap(~feat,labeller = facet_titles)

ggsave("./results_summary/protein_singles_simple_feats_order.png",width = 2500,height=2500,units="px")


ggplot(agg_results %>% dplyr::filter(tax_level == "Order"),aes(x = feat, y = mean_accuracy)) +
  # geom_point(aes(col = Ablation),position = "jitter",size = 2) +
  geom_point(aes(col = Ablation),size = 2,shape = 1) +
  # geom_point(data = noinfo_order,aes(x = feat, y = mean_accuracy),col="grey3",shape = 4) +
  geom_point(data = agg_results_full %>% dplyr::filter(tax_level == "Order"),aes(x = feat,y = max_acc),
             col = "magenta2",shape = 8,size = 1.5)+
  geom_hline(data = noinfo_order,aes(yintercept = mean_accuracy),col="grey3",linetype = "dashed") +
  scale_color_manual(values = palette_OkabeIto_6) + 
  theme_bw() +
  scale_y_continuous(breaks = seq(0,1,0.1))+
  # guides(col = "none") +
  labs(col = "Feature source")+
  theme(strip.background = element_rect(fill="white"),
        axis.text.x = element_text(angle = 315,hjust=0),
        panel.grid.minor = element_blank(),
        legend.position = "top") +
  facet_wrap(~Members,labeller = facet_titles)

ggsave("./results_summary/protein_singles_simple_feats_order_bysubgraph.png",width = 2500,height=2500,units="px")

ggplot(agg_results %>% dplyr::filter(tax_level == "Order" & feat == "dponly"),aes(x = Members, y = mean_accuracy)) +
  # geom_point(aes(col = Ablation),position = "jitter",size = 2,shape = 1) +
  geom_point(aes(col = Ablation,fill=Ablation),size = 3,alpha=0.8,shape = 21) +
  # geom_point(data = noinfo_order,aes(x = feat, y = mean_accuracy),col="grey3",shape = 4) +
  geom_point(data = agg_results_full %>% dplyr::filter(tax_level == "Order" & feat == "dponly"),aes(x = Members,y = max_acc),
             col = "magenta2",shape = 8,size = 1.5)+
  # geom_point(data = noinfo_order %>% dplyr::filter(feat == "dponly"),aes(x = Members, y = mean_accuracy),col="grey3",shape = 4) +
  scale_color_manual(values = palette_OkabeIto_5) + 
  scale_fill_manual(values = palette_OkabeIto_5) + 
  theme_bw() +
  xlab("")+
  ylab("Accuracy")+
  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(0,1))+
  # guides(col = "none") +
  labs(col = "Feature source",fill = "Feature source")+
  theme(strip.background = element_rect(fill="white"),
        axis.text.x = element_text(angle = 315,hjust=0),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        plot.margin = margin(0,1,0,0.2, "cm"))

ggsave("./results_summary/protein_singles_dponly_order.pdf",width = 2500,height=2000,units="px")

ggplot(agg_results %>% dplyr::filter(tax_level == "Order"),aes(x = Members, y = mean_accuracy)) +
  # geom_point(aes(col = Ablation),position = "jitter",size = 2,shape = 1) +
  geom_point(aes(col = Ablation,fill=Ablation),size = 3,alpha=0.4,shape = 21) +
  # geom_point(data = noinfo_order,aes(x = feat, y = mean_accuracy),col="grey3",shape = 4) +
  geom_point(data = agg_results_full %>% dplyr::filter(tax_level == "Order"),aes(x = Members,y = max_acc),
             col = "magenta2",shape = 8,size = 1.5)+
  # geom_point(data = noinfo_order %>% dplyr::filter(feat == "dponly"),aes(x = Members, y = mean_accuracy),col="grey3",shape = 4) +
  scale_color_manual(values = palette_OkabeIto_6) + 
  scale_fill_manual(values = palette_OkabeIto_6) + 
  theme_bw() +
  xlab("")+
  ylab("Accuracy")+  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(0,1))+
  # guides(col = "none") +
  labs(col = "Feature source",fill = "Feature source")+
  theme(strip.background = element_rect(fill="white"),
        axis.text.x = element_text(angle = 315,hjust=0),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        plot.margin = margin(0,2.5,0,0.2, "cm")) +
  facet_wrap(~feat,labeller = facet_titles)

ggsave("./results_summary/protein_singles_simple_feats_order.png",width = 3000,height=2500,units="px")

# DP only 

dp_results <- agg_results %>% 
  filter(tax_level == "Order" & feat == "dponly" & Members != "Reptile viruses" & Members != "Paraavulavirus wisconsinense")

dp_results_full <- agg_results_full %>% 
  filter(tax_level == "Order" & feat == "dponly" & Members != "Reptile viruses" & Members != "Paraavulavirus wisconsinense")

kelly_max_contrast <- c("#FFB300","#A6BDD7","#C10020","#F6768E","#00538A")

palette_check(kelly_max_contrast, plot = T)

ggplot(dp_results %>% dplyr::filter(tax_level == "Order"),aes(x = Members, y = mean_accuracy)) +
  geom_point(aes(col = Ablation,fill=Ablation),size = 3,shape = 23,position = position_dodge(width = 0.4)) +
  geom_point(data = dp_results_full %>% dplyr::filter(tax_level == "Order"),aes(x = Members,y = max_acc),
             col = "black",shape = 8,size = 3)+
  scale_color_manual(values = kelly_max_contrast) + 
  scale_fill_manual(values = kelly_max_contrast) + 
  theme_bw() +
  xlab("")+
  ylab("Accuracy")+  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(0,1))+
  scale_y_continuous(breaks = seq(0.5,1,0.1),limits = c(0.5,1))+
  # guides(col = "none") +
  labs(col = "Feature source",fill = "Feature source")+
  theme(strip.background = element_rect(fill="white"),
        axis.text.x = element_text(angle = 315,hjust=0),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        plot.margin = margin(0,2.5,0,0.2, "cm"))

ggsave("./results_summary/protein_singles_dponly.png",width = 4000,height=4000,units="px",dpi=600)
