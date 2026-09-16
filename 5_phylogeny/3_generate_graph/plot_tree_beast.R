library(ape)
library(tidyverse)
library(RColorBrewer)
library(magrittr)
library(rcartocolor)
library(ggtree)
p_phylo <- ape::read.nexus("../2_beast_phylo/GTR_HMC_HPSTR_tree_50.nxs")

hosts <- read_csv("../../3_host_assignment/all_cluster_hosts.csv")

phylo_tips <- tibble(ref_accessions = p_phylo$tip.label)

tip_hosts <- left_join(phylo_tips, hosts)

tip_hosts %<>%
  mutate(Host_rank = case_when(is.na(Host_rank) == T ~ "Undetermined",
                               .default = Host_rank))

p_phylo$tip.label <- factor(tip_hosts$Host_rank,levels=c("Primates","Rodentia","Carnivora","Chiroptera","Artiodactyla","Perissodactyla","Scandentia","Pholidota","Eulipotyphla",
                                                            "Aves","Reptilia","Fish","Undetermined"))

p_phylo <- ape::drop.tip(p_phylo,"Undetermined")

tip_hosts %<>%
  filter(Host_rank != "Undetermined")

axis_fun <- function() {
  col <- "black"
  for (i in 1:2)
    axis(i, col = col, col.ticks = col, col.axis = col, las = 1)
  box(lty = "19")
}

# set up a colour palette for plotting
palette <- c(carto_pal(12, "Safe"),"grey1")
mycol = palette[p_phylo$tip.label]


p_phylo$tip.label <- tip_hosts$ref_accessions

# p_phylo$tip.label <- as.character(p_phylo$tip.label)

pdf(file="./graph_phylo_plots/host_labels_phylo.pdf",width=16,height=36)
plot(p_phylo, tip.color = mycol); axis_fun()
dev.off()

p_phylo$tip.label <- tip_hosts$ref_accessions

pdf(file="./graph_phylo_plots/accession_labels_phylo.pdf",width=16,height=36)
plot(p_phylo, tip.color = mycol); axis_fun()
dev.off()

ggtree(p_phylo) +
  geom_tiplab(geom = "label",size = 3,color="white",fill = mycol,nudge_x = -0.1) +
  coord_cartesian(clip = "off")

ggsave("./graph_phylo_plots/accession_labels_background_phylo_noundetermined.pdf",width = 16,height=36,units="in")

pdf(file="./graph_phylo_plots/circle_labels_phylo.pdf",width=16,height=36)
plot(p_phylo,show.tip.label=F); axis_fun()
tiplabels(text = NA,pch = 21, bg = mycol, col = mycol,cex=1.8)
dev.off()
