library(tidyverse)
library(writexl)
library(readxl)

read_cat_results <- function(x){
  setwd(x)
  var1 <- as.numeric(str_extract(x,"^\\d+.\\d+"))
  var2 <- as.numeric(str_extract(x,"(?<=_)\\d+.\\d+(?=_)"))
  
  
  temp <- read_csv(paste0("min_seq_id_",var1,"_coverage_",var2,"_stats.csv"))
  temp$params <- paste(var1,var2,sep="_")
  setwd("..")
  return(temp)
}

filelist <- list.files(pattern = "*_clusters$")

all_results <- do.call(rbind,lapply(filelist,read_cat_results))

read_size_metrics <- function(x){
  setwd(x)
  var1 <- as.numeric(str_extract(x,"^\\d+.\\d+"))
  var2 <- as.numeric(str_extract(x,"(?<=_)\\d+.\\d+(?=_)"))
  
  
  temp <- read_csv(paste0("min_seq_id_",var1,"_coverage_",var2,"_cluster_size_metrics.csv"))
  temp$params <- paste(var1,var2,sep="_")
  setwd("..")
  return(temp)
}

all_size_metrics <- do.call(rbind,lapply(filelist,read_size_metrics))

read_size <- function(x){
  setwd(x)
  var1 <- as.numeric(str_extract(x,"^\\d+.\\d+"))
  var2 <- as.numeric(str_extract(x,"(?<=_)\\d+.\\d+(?=_)"))
  
  
  temp <- read_csv(paste0("min_seq_id_",var1,"_coverage_",var2,"_cluster_size.csv"))
  temp$params <- paste(var1,var2,sep="_")
  setwd("..")
  return(temp)
}

all_size <- do.call(rbind,lapply(filelist,read_size))

read_cat_species_raw <- function(x){
  setwd(x)
  var1 <- as.numeric(str_extract(x,"^\\d+.\\d+"))
  var2 <- as.numeric(str_extract(x,"(?<=_)\\d+.\\d+(?=_)"))
  
  
  temp <- read_csv(paste0("min_seq_id_",var1,"_coverage_",var2,"_virus_species_raw.csv"))
  temp$params <- paste(var1,var2,sep="_")
  setwd("..")
  return(temp)
}

raw_species_results <- do.call(rbind,lapply(filelist,read_cat_species_raw))

read_cat_species <- function(x){
  setwd(x)
  var1 <- as.numeric(str_extract(x,"^\\d+.\\d+"))
  var2 <- as.numeric(str_extract(x,"(?<=_)\\d+.\\d+(?=_)"))
  
  
  temp <- read_csv(paste0("min_seq_id_",var1,"_coverage_",var2,"_virus_species_cluster_count.csv"))
  temp$params <- paste(var1,var2,sep="_")
  setwd("..")
  return(temp)
}

all_species_results <- do.call(rbind,lapply(filelist,read_cat_species))

read_cat_centroid_raw <- function(x){
  setwd(x)
  var1 <- as.numeric(str_extract(x,"^\\d+.\\d+"))
  var2 <- as.numeric(str_extract(x,"(?<=_)\\d+.\\d+(?=_)"))
  
  
  temp <- read_csv(paste0("min_seq_id_",var1,"_coverage_",var2,"_centroid_species_raw.csv"))
  temp$params <- paste(var1,var2,sep="_")
  setwd("..")
  return(temp)
}

raw_centroid_results <- do.call(rbind,lapply(filelist,read_cat_centroid_raw))

read_cat_centroid <- function(x){
  setwd(x)
  var1 <- as.numeric(str_extract(x,"^\\d+.\\d+"))
  var2 <- as.numeric(str_extract(x,"(?<=_)\\d+.\\d+(?=_)"))
  
  
  temp <- read_csv(paste0("min_seq_id_",var1,"_coverage_",var2,"_centroid_species_cluster_count.csv"))
  temp$params <- paste(var1,var2,sep="_")
  setwd("..")
  return(temp)
}

all_centroid_results <- do.call(rbind,lapply(filelist,read_cat_centroid))

all_size_metrics$percentage_singletons <- 1-all_size_metrics$percentage_singletons

### NOTE ###
# we add 1 to all counts to allow for plotting on a log y scale (otherwise we have numerous groups with only 1 obs that come out at 0 following log transformation)
# this obviously distorts the relative size of columns so we should be VERY CLEAR ABOUT THIS
# but for this purpose it is useful
# it allows us to assess the chaning shape of the distribution and movement of individual cluster sizes

all_species_density <- tibble(species = NULL,params = NULL)

for(i in 1:nrow(all_species_results)){
  
  all_species_density <- rbind(all_species_density,tibble(species = rep(all_species_results$species[i],times = all_species_results$clusters[i]+1),params = all_species_results$params[i]))

}

all_centroid_density <- tibble(centroids = NULL,params = NULL)

for(i in 1:nrow(all_centroid_results)){
  
  all_centroid_density <- rbind(all_centroid_density,tibble(centroids = rep(all_centroid_results$number_of_centroids[i],times = all_centroid_results$count[i]+1),params = all_centroid_results$params[i]))
  
}


ggplot(all_species_density,aes(x = species)) +
  stat_density(bw=0.05) +
  theme_bw(base_size=12) +
  ylab("Density")+
  xlab("# of viral species represented per cluster")+
  scale_x_log10()+
  facet_wrap(~params)
 
ggsave("species_by_cluster_density.png")

ggplot(all_species_density,aes(species)) +
  geom_histogram() +
  theme_bw(base_size=12) +
  ylab("Count")+
  xlab("# of viral species represented per cluster")+
  scale_y_log10()+
  scale_x_continuous(breaks = 1:8)+
  facet_wrap(~params)

ggsave("species_by_cluster_hist.png")


ggplot(all_centroid_density,aes(x = centroids)) +
  stat_density(bw=0.05) +
  theme_bw(base_size=12) +
  ylab("Density")+
  xlab("# of clusters per centroid species") +
  scale_x_log10()+
  facet_wrap(~params)

ggsave("centroid_species_density.png")


ggplot(all_centroid_density,aes(centroids)) +
  geom_histogram() +
  theme_bw(base_size=12) +
  ylab("Count")+
  xlab("# of clusters per centroid species")+
  scale_y_log10()+
  scale_x_continuous(breaks = c(0,10,seq(from=20,to=80,by=20)))+
  facet_wrap(~params)

ggsave("centroid_species_hist.png")

write_xlsx(all_results,"clustering_tune_full_results.xlsx")

percent_singletons_barplot <- ggplot(all_size_metrics, aes(x = params, y = percentage_singletons, fill = params)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d() +
  labs(title = "Clustering Efficiency", 
       x = "Parameters", y = "% singletons") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

percent_singletons_barplot

ggsave("percent_singletons.png")


total_clusters_barplot <- ggplot(all_size_metrics, aes(x = params, y = total_clusters, fill = params)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d() +
  labs(title = "Total number of clusters", 
       x = "Parameters", y = "# clusters") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

total_clusters_barplot

ggsave("total_clusters.png")


percent_multi_order <- ggplot(all_results, aes(x = params, y = multi_order_prop, fill = params)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d() +
  labs(title = "Proportion of clusters containing viruses with hosts from multiple taxonomic orders", 
       x = "Parameters", y = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

percent_multi_order

ggsave("multi_order.png")


percent_multi_class <- ggplot(all_results, aes(x = params, y = multi_class_prop, fill = params)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d() +
  labs(title = "Proportion of clusters containing viruses with hosts from multiple taxonomic classes", 
       x = "Parameters", y = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

percent_multi_class

ggsave("multi_class.png")


cluster_size_boxplot <- ggplot(all_size, aes(x = params, y = cluster_size, fill = params)) +
  geom_boxplot() +
  scale_y_log10() +
  scale_fill_viridis_d() +
  labs(title = "Cluster Sizes by Parameters",
       x = "Parameters", y = "Cluster Size (log10)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cluster_size_boxplot

ggsave("cluster_size_boxplot.png")


multi_clusts <- raw_centroid_results %>% 
  filter(ref_organism == "Measles morbillivirus" | 
           ref_organism == "Avian orthoavulavirus 1" | 
           ref_organism == "Morbillivirus canis")

write_xlsx(multi_clusts, "big_multi_clusters.xlsx")
