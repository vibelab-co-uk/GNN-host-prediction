library(tidyverse)
library(readxl)
library(viridis)
library(magrittr)

`%nin%` <- Negate(`%in%`)

# read in subgraph results
agg_results <- read_xlsx(paste0("./results_summary/subgraph_agg_performance_best_tune.xlsx"))
roc_results <- read_xlsx(paste0("./results_summary/subgraph_roc_performance_best_tune.xlsx"))
noinfo_results <- read_xlsx(paste0("./results_summary/subgraph_noinfo_performance_best_tune.xlsx"))
agg_results_0 <- read_xlsx(paste0("../../blocked_cv/simple_feats/results_summary/subgraph_agg_performance_best_tune_testonly.xlsx"))

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
  group_by(Block,tax_level,N_samples) %>% 
  summarise(mean_macro = mean(mean_macro,na.rm=T),mean_accuracy = mean(mean_accuracy),sd_accuracy = mean(sd_accuracy))

results_rand$feat <- "rand"

agg_results <- rbind(agg_results,results_rand)

# process random repeats for no reintros
results_0_rand <- agg_results_0 %>% 
  filter(feat %in% c("rand1","rand2","rand3","rand4","rand5"))

agg_results_0 <- agg_results_0 %>% 
  filter(feat %nin% c("rand1","rand2","rand3","rand4","rand5"))

results_0_rand %<>% 
  group_by(Subgraph,tax_level) %>% 
  summarise(mean_macro = mean(mean_macro,na.rm=T),mean_accuracy = mean(mean_accuracy))
  
results_0_rand$feat <- "rand"

agg_results_0 <- rbind(agg_results_0,results_0_rand)

agg_results$feat <- factor(agg_results$feat,levels = c("gconly","dintonly","kmer3only","aaconly","dponly","rand"))
roc_results$feat <- factor(roc_results$feat,levels = c("gconly","dintonly","kmer3only","aaconly","dponly","rand"))
noinfo_results$feat <- factor(noinfo_results$feat,levels = c("gconly","dintonly","kmer3only","aaconly","dponly","rand"))
agg_results_0$feat <- factor(agg_results_0$feat,levels = c("gconly","dintonly","kmer3only","aaconly","dponly","rand"))

agg_results %<>% 
  dplyr::rename("Subgraph" = "Block")

roc_results %<>% 
  dplyr::rename("Subgraph" = "Block")

noinfo_results %<>% 
  dplyr::rename("Subgraph" = "Block")

subgraph_id <- read_xlsx("../../../../5_phylogeny/3_generate_graph/subgraph_virus_classes.xlsx")

best_tune <- read_csv("../../kfold_cv_early_stopping/simple_feats/results_summary/best_tune_acc.csv")

agg_results <- left_join(agg_results,subgraph_id)
roc_results <- left_join(roc_results,subgraph_id)
noinfo_results <- left_join(noinfo_results,subgraph_id)
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

roc_results$Members <- factor(roc_results$Members,levels = subgraph_levels)

noinfo_results$Members <- factor(noinfo_results$Members,levels = subgraph_levels)

# split no info tables into order and class
noinfo_order <- noinfo_results %>% 
  dplyr::filter(tax_level == "Order") %>% 
  distinct()

noinfo_class <- noinfo_results %>% 
  dplyr::filter(tax_level == "Class") %>% 
  distinct()


palette_OkabeIto <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                      "#0072B2", "#D55E00", "#CC79A7", "#999999")


palette_OkabeIto_6 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                        "#CC79A7", "#999999")

agg_results %<>%
  mutate(feat = case_when(feat == "dintonly" ~ "Dinucleotide composition",
                            feat == "aaconly" ~ "Amino acid composition",
                            feat == "dponly" ~ "Dipeptide composition",
                            feat == "kmer3only" ~ "3mer composition",
                            feat == "gconly" ~ "GC %",
                            feat == "rand" ~ "Random pseudofeature"))

agg_results$feat <- factor(agg_results$feat,levels = c("GC %","Dinucleotide composition","3mer composition",
                                                       "Amino acid composition","Dipeptide composition","Random pseudofeature"))

agg_results %<>% 
  filter(tax_level == "Order")

ggplot(agg_results,aes(x = N_samples, y = mean_accuracy,col = feat)) +
  # geom_point(size = 2) +
  geom_line()+
  geom_pointrange(data = agg_results,
                  aes(x = N_samples,y=mean_accuracy,ymin = mean_accuracy-sd_accuracy,ymax = mean_accuracy+sd_accuracy),
                  size = 0.5,shape = 4)+
  scale_color_manual(values = palette_OkabeIto_6) +
  theme_bw() +
  xlab("In-block train samples")+
  ylab("Mean accuracy")+
  scale_y_continuous(limits = c(0,1))+
  theme(strip.background = element_rect(fill="white"),
        panel.grid.minor = element_blank())+
        # legend.position = "top") +
  guides(col=guide_legend(title="Feature set")) +
  facet_wrap(~Members,nrow=2)
# facet_wrap(~Subgraph,labeller = facet_titles,nrow=2)

ggsave("./results_summary/reintro_accuracy_gain.png",width = 3500,height=1500,units="px")

# plot with the blast benchmark
palette_OkabeIto_6 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                        "#CC79A7", "#999999","black")

blast_results$feat <- "BLAST benchmark"

agg_results %<>%
  select(mean_accuracy,sd_accuracy,feat,Subgraph,N_samples,Members)

agg_results <- rbind(agg_results,blast_results)

agg_results$feat <- factor(agg_results$feat,levels = c("GC %","Dinucleotide composition","3mer composition",
                                                       "Amino acid composition","Dipeptide composition","Random pseudofeature","BLAST benchmark"))

ggplot(agg_results,aes(x = N_samples, y = mean_accuracy,col = feat)) +
  # geom_point(size = 2) +
  geom_line()+
  geom_pointrange(data = agg_results,
                  aes(x = N_samples,y=mean_accuracy,ymin = mean_accuracy-sd_accuracy,ymax = mean_accuracy+sd_accuracy),
                  size = 0.5,shape = 4)+
  scale_color_manual(values = palette_OkabeIto_6) +
  theme_bw() +
  xlab("In-block train samples")+
  ylab("Mean accuracy")+
  scale_y_continuous(limits = c(0,1))+
  theme(strip.background = element_rect(fill="white"),
        panel.grid.minor = element_blank())+
  # legend.position = "top") +
  guides(col=guide_legend(title="Feature set")) +
  facet_wrap(~Members,nrow=2)

ggsave("./results_summary/reintro_accuracy_gain_wBLAST.png",width = 3500,height=1500,units="px")