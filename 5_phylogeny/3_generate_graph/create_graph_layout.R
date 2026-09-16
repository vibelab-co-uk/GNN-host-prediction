library(igraph)
library(ape)
library(tidyverse)
library(writexl)
library(magrittr)
library(RColorBrewer)
library(rcartocolor)

### READ ME ###
################################################################################

# this script was used to generate a manual graph layout for graph visualisation
# the resulting layout is saved as "final_manual_layout.csv" and this script should not need to be run again

################################################################################


# read in maximum likelihood phylogeny
p_phylo <- ape::read.nexus("../2_beast_phylo/GTR_HMC_HPSTR_tree_50.nxs")

hosts <- read_csv("../../3_host_assignment/all_cluster_hosts.csv")

# create tibble of tip labels and bind in hosts
phylo_tips <- tibble(ref_accessions = p_phylo$tip.label)

tip_hosts <- left_join(phylo_tips, hosts)

tip_hosts %<>%
  mutate(Host_rank = case_when(is.na(Host_rank) == T ~ "Undetermined",
                               .default = Host_rank))

# label the phylogeny with hosts and drop tips with undetermined hosts
p_phylo$tip.label <- factor(tip_hosts$Host_rank,levels=c("Primates","Rodentia","Carnivora","Chiroptera","Artiodactyla","Perissodactyla","Scandentia","Pholidota","Eulipotyphla",
                                                         "Aves","Reptilia","Fish","Undetermined"))

p_phylo <- ape::drop.tip(p_phylo,"Undetermined")

tip_hosts %<>%
  filter(Host_rank != "Undetermined")

pat_distances <- ape::cophenetic.phylo(p_phylo)

# generate a graph from this data -----------------------------------------
#_______________________________________________________________________________
#_______________________________________________________________________________

# it already functions as the adjacency matrix of a fully connected graph
g <- graph_from_adjacency_matrix(pat_distances, mode = "max", weighted = T)

# set up a colour palette for plotting
# palette <- c(carto_pal(12, "Safe"),"grey1")
palette <- carto_pal(12, "Safe")

mycol = palette[p_phylo$tip.label]

# generate an mst graph ---------------------------------------------------
mg <- igraph::mst(g, weights = igraph::E(g)$weight) # generate minimum spanning tree

# write_csv(as.data.frame(layout_kk),"./graph_phylo_plots/kk_layout1_logscaled_arbitrary_norm.csv")
layout_kk <- read_csv("./graph_phylo_plots/kk_layout1_logscaled_arbitrary_norm.csv")
layout_tbl <- as_tibble(layout_kk)

colnames(layout_tbl) <- c("x_coord","y_coord")

layout_tbl$host <- p_phylo$tip.label
layout_tbl$accession <- tip_hosts$ref_accessions

# we now manually fix overlaps while maintaining edge lengths
# left overlapping cluster pivots from PV102958.1
# we will pivot this to the right while maintaining edge length and then move all downstream nodes to maintain position relative to this node
left_cluster_ref <- "PV102958.1"
left_cluster_members <- c("OL441592.1","MK037461.1","KC802221.1","MW557651.1",
                          "OP628167.1","MZ328285.1","JQ697837.1","MK188747.1",
                          "MK088517.1","AB924122.1","MN604235.1")

bottom_cluster_ref <- "OM030315.1"
# bottom_cluster_members <- c("OM030316.1","PQ140950.1","OK623354.1","OK623353.1")
bottom_cluster_members <- c("OK623354.1","OK623353.1")

right_cluster_incidental_movers <- c("EU782025.1","EU403085.1")

right_cluster_ref <- "MK677433.1"
right_cluster_members <- c("MZ802808.1","FJ231524.1","PX093055.1","MF448515.1",
                           "HM755886.2","JQ886184.1","MH844488.1","MH844489.1",
                           "PQ300082.1","MW338847.1","MZ852789.1","GU206351.1",
                           "MZ351191.1","PV643989.1")

# right_cluster_top_branch <- c("MF448515.1","HM755886.2","JQ886184.1","MH844488.1",
#                               "MH844489.1","PQ300082.1")

# right_cluster_bottom_branch <- c("MZ852789.1","GU206351.1","MZ351191.1","PV643989.1")


### fix bottom cluster------------------------------------------------
# this cluster can just be moved in the x co-ordinate, y co-ordinates can remain the same
# flip pivot to other side of its neighbour, maintaining edge length
pivot <- layout_tbl %>% filter(accession == "OR713879.1") %>% select(x_coord,y_coord)
ref <- layout_tbl %>% filter(accession == bottom_cluster_ref) %>% select(x_coord,y_coord)
members <- layout_tbl %>% filter(accession %in% bottom_cluster_members)

# calculate the x distance of all members from the ref's original position
members %<>%
  mutate(old_dist = x_coord - ref$x_coord)

# calculate distance to move the ref relative to the pivot
ref_dist <- 2*(pivot$x_coord - ref$x_coord)

# move the ref relative to the pivot
layout_tbl %<>%
  mutate(x_coord = case_when(accession == bottom_cluster_ref ~ x_coord + ref_dist,
                             .default = x_coord))

# we move all members by the amount the ref moved + double their original distance from the ref
for(i in 1:nrow(members)){
  layout_tbl %<>%
    mutate(x_coord = case_when(accession == members$accession[i] ~ x_coord + (ref_dist + abs(2*members$old_dist[i])),
                               .default = x_coord))
}

### fix left cluster------------------------------------------------
# this cluster can just be moved in the x co-ordinate, y co-ordinates can remain the same
# flip ref to other side of its neighbour (the pivot), maintaining edge length
pivot <- layout_tbl %>% filter(accession == "MT511667.1") %>% select(x_coord,y_coord)
ref <- layout_tbl %>% filter(accession == left_cluster_ref) %>% select(x_coord,y_coord)
members <- layout_tbl %>% filter(accession %in% left_cluster_members)

# calculate the distance of all members from the ref's original position
members %<>%
  mutate(old_dist = x_coord - ref$x_coord)

# calculate distance to move the ref relative to the pivot
ref_dist <- 2*(pivot$x_coord - ref$x_coord)

# move the ref relative to the pivot
layout_tbl %<>%
  mutate(x_coord = case_when(accession == left_cluster_ref ~ x_coord + ref_dist,
                             .default = x_coord))

# we move all members by the amount the ref moved + double their original distance from the ref
for(i in 1:nrow(members)){
  layout_tbl %<>%
    mutate(x_coord = case_when(accession == members$accession[i] & members$old_dist[i] <0 ~ x_coord + (ref_dist + abs(2*members$old_dist[i])),
                               accession == members$accession[i] & members$old_dist[i] >= 0 ~ x_coord + ref_dist,
                               .default = x_coord))
}

# move the ref again to make a bit more space
# set the x co-ordinate position we want, then calculate the y-coordinate needed to maintain distance

# calculate current distance between ref and pivot
# now find the missing value, we have both co-ordinates of one position, the distance and one co-ordinate of the other
# we can use pythagoras to find the value we want

ref_pivot_dist <- sqrt((ref$x_coord - pivot$x_coord)^2 +
                         (ref$y_coord - pivot$y_coord)^2)

new_x <- pivot$x_coord # we'll move to directly above the pivot

# 6 calculate the new y coord
new_y = sqrt(ref_pivot_dist^2 - (new_x - pivot$x_coord)^2) + pivot$y_coord 

layout_tbl %<>%
  mutate(x_coord = case_when(accession == left_cluster_ref ~ new_x,
                             .default = x_coord),
         y_coord = case_when(accession == left_cluster_ref ~ new_y,
                             .default = y_coord))

new_ref <- layout_tbl %>% filter(accession == left_cluster_ref) %>% select(x_coord,y_coord)

move <- ref - new_ref

ref_pivot_dist_check <- sqrt((new_ref$x_coord - pivot$x_coord)^2 +
                         (new_ref$y_coord - pivot$y_coord)^2)

ref_pivot_dist == ref_pivot_dist_check 

# then move all members to maintain relative position to the ref
for(i in 1:nrow(members)){
  layout_tbl %<>%
    mutate(x_coord = case_when(accession == members$accession[i] ~ x_coord + move$x_coord,
                               .default = x_coord),
           y_coord = case_when(accession == members$accession[i] ~ y_coord - move$y_coord,
                               .default = y_coord))
}

# we need a secondary pivot to avoid overlaps
left_cluster_ref2 <- "MK188747.1"
left_cluster_members2 <- c("MK088517.1","AB924122.1","MN604235.1")

pivot <- layout_tbl %>% filter(accession == "OP628167.1") %>% select(x_coord,y_coord)
ref <- layout_tbl %>% filter(accession == left_cluster_ref2) %>% select(x_coord,y_coord)
members <- layout_tbl %>% filter(accession %in% left_cluster_members2)

# set the x co-ordinate position we want, then calculate the y-coordinate needed to maintain distance

# calculate current distance between ref and pivot
# now find the missing value, we have both co-ordinates of one position, the distance and one co-ordinate of the other
# we can use pythagoras to find the value we want
ref_pivot_dist <- sqrt((ref$x_coord - pivot$x_coord)^2 +
                         (ref$y_coord - pivot$y_coord)^2)

new_x <- -2

# solve new y coordinate of ref maintaining this distance
new_y = sqrt(abs(ref_pivot_dist^2 - (new_x - pivot$x_coord)^2)) + pivot$y_coord 

layout_tbl %<>%
  mutate(x_coord = case_when(accession == left_cluster_ref2 ~ new_x,
                             .default = x_coord),
         y_coord = case_when(accession == left_cluster_ref2 ~ new_y,
                             .default = y_coord))

new_ref <- layout_tbl %>% filter(accession == left_cluster_ref2) %>% select(x_coord,y_coord)

move <- ref - new_ref

# then move all members to maintain relative position to the ref
for(i in 1:nrow(members)){
  layout_tbl %<>%
    mutate(x_coord = case_when(accession == members$accession[i] ~ x_coord - move$x_coord,
                               .default = x_coord),
           y_coord = case_when(accession == members$accession[i] ~ y_coord - move$y_coord,
                               .default = y_coord))
}

# fix right cluster -------------------------------------------------------
# first move a couple out the way
pivot <- layout_tbl %>% filter(accession == "MH044693.1") %>% select(x_coord,y_coord)
ref <- layout_tbl %>% filter(accession == "EU782025.1") %>% select(x_coord,y_coord)
members <- layout_tbl %>% filter(accession %in% "EU403085.1")

ref_pivot_dist <- sqrt((ref$x_coord - pivot$x_coord)^2 +
                         (ref$y_coord - pivot$y_coord)^2)

new_x <- 2.2

# solve new y coordinate of ref maintaining this distance
new_y = -sqrt(ref_pivot_dist^2 - (new_x - pivot$x_coord)^2) + pivot$y_coord 

layout_tbl %<>%
  mutate(x_coord = case_when(accession == "EU782025.1" ~ new_x,
                             .default = x_coord),
         y_coord = case_when(accession == "EU782025.1" ~ new_y,
                             .default = y_coord))

new_ref <- layout_tbl %>% filter(accession == "EU782025.1") %>% select(x_coord,y_coord)

move <- ref - new_ref

# then move all members to maintain relative position to the ref
for(i in 1:nrow(members)){
  layout_tbl %<>%
    mutate(x_coord = case_when(accession == members$accession[i] ~ x_coord - move$x_coord,
                               .default = x_coord),
           y_coord = case_when(accession == members$accession[i] ~ y_coord - move$y_coord,
                               .default = y_coord))
}

# now pivot the big chain
pivot <- layout_tbl %>% filter(accession == "OR367435.1") %>% select(x_coord,y_coord)
ref <- layout_tbl %>% filter(accession == right_cluster_ref) %>% select(x_coord,y_coord)
members <- layout_tbl %>% filter(accession %in% right_cluster_members)

ref_pivot_dist <- sqrt((ref$x_coord - pivot$x_coord)^2 +
                         (ref$y_coord - pivot$y_coord)^2)

new_y <- 1.2

# solve new y coordinate of ref maintaining this distance
new_x = sqrt(abs(ref_pivot_dist^2 - (new_y - pivot$y_coord)^2)) + pivot$x_coord 

layout_tbl %<>%
  mutate(x_coord = case_when(accession == right_cluster_ref ~ new_x,
                             .default = x_coord),
         y_coord = case_when(accession == right_cluster_ref ~ new_y,
                             .default = y_coord))

new_ref <- layout_tbl %>% filter(accession == right_cluster_ref) %>% select(x_coord,y_coord)

move <- ref - new_ref

# then move all members to maintain relative position to the ref
for(i in 1:nrow(members)){
  layout_tbl %<>%
    mutate(x_coord = case_when(accession == members$accession[i] ~ x_coord - move$x_coord,
                               .default = x_coord),
           y_coord = case_when(accession == members$accession[i] ~ y_coord - move$y_coord,
                               .default = y_coord))
}

# second pivot
right_cluster_ref2 <- "FJ231524.1"
right_cluster_members2 <- c("PX093055.1","MF448515.1",
                            "HM755886.2","JQ886184.1","MH844488.1","MH844489.1",
                            "PQ300082.1","MW338847.1","MZ852789.1","GU206351.1",
                            "MZ351191.1","PV643989.1")
  
pivot <- layout_tbl %>% filter(accession == "MZ802808.1") %>% select(x_coord,y_coord)
ref <- layout_tbl %>% filter(accession == right_cluster_ref2) %>% select(x_coord,y_coord)
members <- layout_tbl %>% filter(accession %in% right_cluster_members2)

ref_pivot_dist <- sqrt((ref$x_coord - pivot$x_coord)^2 +
                         (ref$y_coord - pivot$y_coord)^2)

new_y <- 0.9

# solve new y coordinate of ref maintaining this distance
new_x = sqrt(ref_pivot_dist^2 - (new_y - pivot$y_coord)^2) + pivot$x_coord 

layout_tbl %<>%
  mutate(x_coord = case_when(accession == right_cluster_ref2 ~ new_x,
                             .default = x_coord),
         y_coord = case_when(accession == right_cluster_ref2 ~ new_y,
                             .default = y_coord))

new_ref <- layout_tbl %>% filter(accession == right_cluster_ref2) %>% select(x_coord,y_coord)

move <- ref - new_ref

# then move all members to maintain relative position to the ref
for(i in 1:nrow(members)){
  layout_tbl %<>%
    mutate(x_coord = case_when(accession == members$accession[i] ~ x_coord - move$x_coord,
                               .default = x_coord),
           y_coord = case_when(accession == members$accession[i] ~ y_coord - move$y_coord,
                               .default = y_coord))
}

# third pivot
right_cluster_ref2 <- "PX093055.1"
right_cluster_members2 <- c("MF448515.1",
                            "HM755886.2","JQ886184.1","MH844488.1","MH844489.1",
                            "PQ300082.1","MW338847.1","MZ852789.1","GU206351.1",
                            "MZ351191.1","PV643989.1")

pivot <- layout_tbl %>% filter(accession == "FJ231524.1") %>% select(x_coord,y_coord)
ref <- layout_tbl %>% filter(accession == right_cluster_ref2) %>% select(x_coord,y_coord)
members <- layout_tbl %>% filter(accession %in% right_cluster_members2)

ref_pivot_dist <- sqrt((ref$x_coord - pivot$x_coord)^2 +
                         (ref$y_coord - pivot$y_coord)^2)

new_y <- 0.5

# solve new y coordinate of ref maintaining this distance
new_x = sqrt(ref_pivot_dist^2 - (new_y - pivot$y_coord)^2) + pivot$x_coord 

layout_tbl %<>%
  mutate(x_coord = case_when(accession == right_cluster_ref2 ~ new_x,
                             .default = x_coord),
         y_coord = case_when(accession == right_cluster_ref2 ~ new_y,
                             .default = y_coord))

new_ref <- layout_tbl %>% filter(accession == right_cluster_ref2) %>% select(x_coord,y_coord)

move <- ref - new_ref

# then move all members to maintain relative position to the ref
for(i in 1:nrow(members)){
  layout_tbl %<>%
    mutate(x_coord = case_when(accession == members$accession[i] ~ x_coord - move$x_coord,
                               .default = x_coord),
           y_coord = case_when(accession == members$accession[i] ~ y_coord - move$y_coord,
                               .default = y_coord))
}

# fourth pivot
right_cluster_ref2 <- "MF448515.1"
right_cluster_members2 <- c("HM755886.2","JQ886184.1","MH844488.1","MH844489.1",
                            "PQ300082.1")

pivot <- layout_tbl %>% filter(accession == "PX093055.1") %>% select(x_coord,y_coord)
ref <- layout_tbl %>% filter(accession == right_cluster_ref2) %>% select(x_coord,y_coord)
members <- layout_tbl %>% filter(accession %in% right_cluster_members2)

ref_pivot_dist <- sqrt((ref$x_coord - pivot$x_coord)^2 +
                         (ref$y_coord - pivot$y_coord)^2)

new_y <- 0.6

# solve new y coordinate of ref maintaining this distance
new_x = sqrt(ref_pivot_dist^2 - (new_y - pivot$y_coord)^2) + pivot$x_coord 

layout_tbl %<>%
  mutate(x_coord = case_when(accession == right_cluster_ref2 ~ new_x,
                             .default = x_coord),
         y_coord = case_when(accession == right_cluster_ref2 ~ new_y,
                             .default = y_coord))

new_ref <- layout_tbl %>% filter(accession == right_cluster_ref2) %>% select(x_coord,y_coord)

move <- ref - new_ref

# then move all members to maintain relative position to the ref
for(i in 1:nrow(members)){
  layout_tbl %<>%
    mutate(x_coord = case_when(accession == members$accession[i] ~ x_coord - move$x_coord,
                               .default = x_coord),
           y_coord = case_when(accession == members$accession[i] ~ y_coord - move$y_coord,
                               .default = y_coord))
}

# fifth pivot
right_cluster_ref2 <- "HM755886.2"
right_cluster_members2 <- c("JQ886184.1","MH844488.1","MH844489.1",
                            "PQ300082.1")

pivot <- layout_tbl %>% filter(accession == "MF448515.1") %>% select(x_coord,y_coord)
ref <- layout_tbl %>% filter(accession == right_cluster_ref2) %>% select(x_coord,y_coord)
members <- layout_tbl %>% filter(accession %in% right_cluster_members2)

ref_pivot_dist <- sqrt((ref$x_coord - pivot$x_coord)^2 +
                         (ref$y_coord - pivot$y_coord)^2)

new_y <- 0.7

# solve new y coordinate of ref maintaining this distance
new_x = sqrt(ref_pivot_dist^2 - (new_y - pivot$y_coord)^2) + pivot$x_coord 

layout_tbl %<>%
  mutate(x_coord = case_when(accession == right_cluster_ref2 ~ new_x,
                             .default = x_coord),
         y_coord = case_when(accession == right_cluster_ref2 ~ new_y,
                             .default = y_coord))

new_ref <- layout_tbl %>% filter(accession == right_cluster_ref2) %>% select(x_coord,y_coord)

move <- ref - new_ref

# then move all members to maintain relative position to the ref
for(i in 1:nrow(members)){
  layout_tbl %<>%
    mutate(x_coord = case_when(accession == members$accession[i] ~ x_coord - move$x_coord,
                               .default = x_coord),
           y_coord = case_when(accession == members$accession[i] ~ y_coord - move$y_coord,
                               .default = y_coord))
}

# sixth pivot
right_cluster_ref2 <- "MH844489.1"
right_cluster_members2 <- c("MH844488.1","PQ300082.1")

pivot <- layout_tbl %>% filter(accession == "HM755886.2") %>% select(x_coord,y_coord)
ref <- layout_tbl %>% filter(accession == right_cluster_ref2) %>% select(x_coord,y_coord)
members <- layout_tbl %>% filter(accession %in% right_cluster_members2)

ref_pivot_dist <- sqrt((ref$x_coord - pivot$x_coord)^2 +
                         (ref$y_coord - pivot$y_coord)^2)

new_y <- 0.2

# solve new y coordinate of ref maintaining this distance
new_x = sqrt(ref_pivot_dist^2 - (new_y - pivot$y_coord)^2) + pivot$x_coord 

layout_tbl %<>%
  mutate(x_coord = case_when(accession == right_cluster_ref2 ~ new_x,
                             .default = x_coord),
         y_coord = case_when(accession == right_cluster_ref2 ~ new_y,
                             .default = y_coord))

new_ref <- layout_tbl %>% filter(accession == right_cluster_ref2) %>% select(x_coord,y_coord)

move <- ref - new_ref

# then move all members to maintain relative position to the ref
for(i in 1:nrow(members)){
  layout_tbl %<>%
    mutate(x_coord = case_when(accession == members$accession[i] ~ x_coord - move$x_coord,
                               .default = x_coord),
           y_coord = case_when(accession == members$accession[i] ~ y_coord - move$y_coord,
                               .default = y_coord))
}

# seventh pivot
right_cluster_ref2 <- "MH844488.1"
right_cluster_members2 <- c("PQ300082.1")

pivot <- layout_tbl %>% filter(accession == "MH844489.1") %>% select(x_coord,y_coord)
ref <- layout_tbl %>% filter(accession == right_cluster_ref2) %>% select(x_coord,y_coord)
members <- layout_tbl %>% filter(accession %in% right_cluster_members2)

# calculate the y distance of all members from the ref's original position
members %<>%
  mutate(old_dist = y_coord - ref$y_coord)

# calculate distance to move the ref relative to the pivot
ref_dist <- 2*(pivot$y_coord - ref$y_coord)

# move the ref relative to the pivot
layout_tbl %<>%
  mutate(y_coord = case_when(accession == right_cluster_ref2 ~ y_coord + ref_dist,
                             .default = y_coord))

# we move all members by the amount the ref moved + double their original distance from the ref
for(i in 1:nrow(members)){
  layout_tbl %<>%
    mutate(y_coord = case_when(accession == members$accession[i] ~ y_coord + (ref_dist - abs(2*members$old_dist[i])),
                               .default = y_coord))
}

# one off fixes -----------------------------------------------------------

pivot <- layout_tbl %>% filter(accession == "MZ328277.1") %>% select(x_coord,y_coord)
ref <- layout_tbl %>% filter(accession == "AF079780.2") %>% select(x_coord,y_coord)
members <- layout_tbl %>% filter(accession == "ON861833.1")

ref_pivot_dist <- sqrt((ref$x_coord - pivot$x_coord)^2 +
                         (ref$y_coord - pivot$y_coord)^2)

new_y <- -4.1

# solve new y coordinate of ref maintaining this distance
new_x = sqrt(abs(ref_pivot_dist^2 - (new_y - pivot$y_coord)^2)) + pivot$x_coord 

layout_tbl %<>%
  mutate(x_coord = case_when(accession == "AF079780.2" ~ new_x,
                             .default = x_coord),
         y_coord = case_when(accession == "AF079780.2" ~ new_y,
                             .default = y_coord))

new_ref <- layout_tbl %>% filter(accession == "AF079780.2") %>% select(x_coord,y_coord)

move <- ref - new_ref

# then move all members to maintain relative position to the ref
for(i in 1:nrow(members)){
  layout_tbl %<>%
    mutate(x_coord = case_when(accession == members$accession[i] ~ x_coord - move$x_coord,
                               .default = x_coord),
           y_coord = case_when(accession == members$accession[i] ~ y_coord - move$y_coord,
                               .default = y_coord))
}

# one off fix 2
pivot <- layout_tbl %>% filter(accession == "MH892405.1") %>% select(x_coord,y_coord)
ref <- layout_tbl %>% filter(accession == "MH972568.1") %>% select(x_coord,y_coord)

ref_pivot_dist <- sqrt((ref$x_coord - pivot$x_coord)^2 +
                         (ref$y_coord - pivot$y_coord)^2)

new_y <- 2

# solve new y coordinate of ref maintaining this distance
new_x = sqrt(ref_pivot_dist^2 - (new_y - pivot$y_coord)^2) + pivot$x_coord 

layout_tbl %<>%
  mutate(x_coord = case_when(accession == "MH972568.1" ~ new_x,
                             .default = x_coord),
         y_coord = case_when(accession == "MH972568.1" ~ new_y,
                             .default = y_coord))

layout_kk <- layout_tbl %>% select(x_coord,y_coord)
layout_kk <- as.matrix(layout_kk)

# save final layout and plot ----------------------------------------------
write_csv(layout_tbl,"./graph_phylo_plots/final_manual_layout.csv")

############
############

png(filename="./graph_phylo_plots/kk_manual_host_graph_wlegend.png",width=2000,heigh=4000,units="px")
plot(mg, vertex.size=3, vertex.label=NA,vertex.color=mycol,layout = layout_kk,edge.width = 3)
legend('topright',legend=levels(factor(tip_hosts$Host_rank,levels=c("Primates","Rodentia","Carnivora","Chiroptera","Artiodactyla","Perissodactyla","Scandentia","Pholidota","Eulipotyphla",
                                                          "Aves","Reptilia","Fish"))),col = palette,pch=16,cex=3)
dev.off()


png(filename="./graph_phylo_plots/kk_manual_host_graph.png",width=2000,heigh=4000,units="px")
plot(mg, vertex.size=3, vertex.label=NA,vertex.color=mycol,layout = layout_kk,edge.width = 3)
dev.off()

png(filename="./graph_phylo_plots/kk_manual_host_graph_labelled.png",width=2000,heigh=4000,units="px")
plot(mg, vertex.size=3, vertex.color=mycol,layout = layout_kk,edge.width = 3)
dev.off()

# relabel tips with accessions instead of hosts
png(filename="./graph_phylo_plots/kk_manual_host_graph_accessionlabelled.png",width=2000,heigh=4000,units="px")
plot(mg, vertex.size=3, vertex.label=tip_hosts$ref_accessions,vertex.color=mycol,layout = layout_kk,edge.width = 3)
dev.off()
