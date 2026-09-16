library(Biostrings)

# short script to rename MSAs with only accessions so they are consistent for BEAST
N_orfs <- readDNAStringSet("../1_msa/cluster_centroids_aligned_ginsi_nucleocapsid_orfs.fasta")
M_orfs <- readDNAStringSet("../1_msa/cluster_centroids_aligned_ginsi_matrix_orfs.fasta")
F_orfs <- readDNAStringSet("../1_msa/cluster_centroids_aligned_ginsi_fusion_orfs.fasta")
HN_orfs <- readDNAStringSet("../1_msa/cluster_centroids_aligned_ginsi_attachment_orfs.fasta")
L_orfs <- readDNAStringSet("../1_msa/cluster_centroids_aligned_ginsi_polymerase_orfs.fasta")

names(N_orfs) <- str_split_i(names(N_orfs),pattern = "\\|",i=1) # rename tips with only sequence IDs
names(M_orfs) <- str_split_i(names(M_orfs),pattern = "\\|",i=1) # rename tips with only sequence IDs
names(F_orfs) <- str_split_i(names(F_orfs),pattern = "\\|",i=1) # rename tips with only sequence IDs
names(HN_orfs) <- str_split_i(names(HN_orfs),pattern = "\\|",i=1) # rename tips with only sequence IDs
names(L_orfs) <- str_split_i(names(L_orfs),pattern = "\\|",i=1) # rename tips with only sequence IDs

N_names <- names(N_orfs)
M_names <- names(M_orfs)
F_names <- names(F_orfs)
HN_names <- names(HN_orfs)
L_names <- names(L_orfs)

writeXStringSet(N_orfs,"N_orfs_aligned.fasta")
writeXStringSet(M_orfs,"M_orfs_aligned.fasta")
writeXStringSet(F_orfs,"F_orfs_aligned.fasta")
writeXStringSet(HN_orfs,"HN_orfs_aligned.fasta")
writeXStringSet(L_orfs,"L_orfs_aligned.fasta")
