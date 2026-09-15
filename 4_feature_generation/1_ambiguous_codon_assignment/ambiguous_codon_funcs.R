# script containing functions required for the 'resolve_ambiguous_codons.R' script
### define functions used to calculate consensus codons for positions in centroid sequences with ambiguous amino acids ###

consensus_finder <- function(query,target,sites){
  
  target_orfs <- orf_table %>% 
    dplyr::filter(seqnames == str_split_i(names(target),"\\s\\|",i = 1) & Protein %in% sites$Protein)
  
  target_orfs$ref_id <- paste(target_orfs$seqnames,target_orfs$Protein,sep="|")
  
  # calculate the start position of the ambiguous codons in the whole nucleotide sequence
  sites %<>%
    mutate(nuc_start = nuc_start + target_orfs$start, nuc_end = nuc_end + target_orfs$start)
  
  alignment <- pwalign::pairwiseAlignment(target,query)
  
  alignment <- c(pwalign::alignedPattern(alignment),pwalign::alignedSubject(alignment))
  
  ambig_align <- subseq(alignment,start=sites$nuc_start,end=sites$nuc_end)[2][[1]]
  
  return(ambig_align)
}

align_ambig <- function(target_site,sites,target_seqs,target_fasta){
  
  target_site <- sites[target_site,]
  
  member_fasta <- ambig_fasta[grepl(paste(target_seqs$members,collapse = "|"),names(ambig_fasta))]
  
  member_fasta <- head(member_fasta)
  
  aligned_codons <- lapply(member_fasta,consensus_finder,target=target_fasta,sites=target_site)
  
  aligned_codons <- DNAStringSet(aligned_codons)
  
  return(aligned_codons) 
}

resolve_ambig <- function(x){
  
  target_seqs <- clusters %>% 
    dplyr::filter(ref_accessions == x) %>% 
    dplyr::filter(members != x)
  
  # if this is a singleton cluster, we have to remove the ambiguous amino acid
  if(nrow(target_seqs) == 0){
    message(paste0("Sequence ",x," processed; no other cluster sequences available"))
    return(NULL)
  } else {
    target_fasta <- ambig_fasta[grepl(x,names(ambig_fasta))]
    
    # filter for the target ambiguous sites and orfs
    target_ambig_sites <- ambig_sites %>% 
      separate_wider_delim(ref_id,delim="|",names=c("seqnames","Protein"),cols_remove = F) %>% 
      dplyr::filter(seqnames == x)
    
    # align each member of the cluster to the centroid sequence and extract the codon at the relevant ambiguous positions
    aligned_codons <- lapply(1:nrow(target_ambig_sites),align_ambig,sites=target_ambig_sites,target_seqs=target_seqs,target_fasta=target_fasta)
    
    for(i in 1:length(aligned_codons)){
      consensus <- consensusMatrix(aligned_codons[[i]])
      
      consensus_codon <- DNAString(paste0(names(which.max(consensus[,1])),names(which.max(consensus[,2])),names(which.max(consensus[,3]))))
      
      target_ambig_sites$consensus_codon[i] <- as.character(consensus_codon)
      target_ambig_sites$consensus_aa[i] <- as.character(translate(consensus_codon,if.fuzzy.codon = "solve"))
    }
    
    message(paste0("Sequence ",x," processed succesfully"))
    
    return(target_ambig_sites)
  }
}
