library(ape)
library(tidyverse)
library(RColorBrewer)
library(magrittr)
library(rcartocolor)

p_phylo <- ape::read.nexus("../2_beast_phylo/L_phylogeny/L_orfs_GTR_Strict_HIPSTR_summary.nxs")

hosts <- read_csv("../../3_host_assignment/all_cluster_hosts.csv")

phylo_tips <- tibble(ref_accessions = p_phylo$tip.label)

tip_hosts <- left_join(phylo_tips, hosts)

tip_hosts %<>%
  mutate(Host_rank = case_when(is.na(Host_rank) == T ~ "Undetermined",
                               .default = Host_rank))

p_phylo$tip.label <- factor(tip_hosts$Host_rank,levels=c("Primates","Rodentia","Carnivora","Chiroptera","Artiodactyla","Perissodactyla","Scandentia","Pholidota","Eulipotyphla",
                                                            "Aves","Reptilia","Fish","Undetermined"))

axis_fun <- function() {
  col <- "black"
  for (i in 1:2)
    axis(i, col = col, col.ticks = col, col.axis = col, las = 1)
  box(lty = "19")
}

# set up a colour palette for plotting
palette <- c(carto_pal(12, "Safe"),"grey1")
mycol = palette[p_phylo$tip.label]

p_phylo$tip.label <- as.character(p_phylo$tip.label)

pdf(file="./graph_phylo_plots/host_labels_phylo_Lonly.pdf",width=16,height=36)
plot(p_phylo, tip.color = mycol); axis_fun()
dev.off()

p_phylo$tip.label <- tip_hosts$ref_accessions

pdf(file="./graph_phylo_plots/accession_labels_phylo_Lonly.pdf",width=16,height=36)
plot(p_phylo, tip.color = mycol); axis_fun()
dev.off()
