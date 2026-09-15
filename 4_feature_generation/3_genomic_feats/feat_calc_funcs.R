# script containing functions required for the 'process_calc_feats_orfs.R' script
### define a function to apply functions to every combination of multiple arguments ###

cmapply <- function(FUN, ..., MoreArgs = NULL){  
  # expand a grid of all argument combinations  
  l <- expand.grid(..., stringsAsFactors=FALSE)    
  # apply the function  
  .mapply(FUN=FUN, dots=unname(l), MoreArgs = MoreArgs)}

### GENOME FEATS ###

calc_composition_prop_gc <- function(x, overlap = TRUE){
  
  df <- cbind(data.frame(
    nuc_id = x %>% names(),
    x %>% letterFrequency("GC", as.prob = TRUE) * 100 %>% as.vector())) %>% 
    rename_at(vars(G.C), ~"GC_content")# Calculate % GC content
  return(df)
  
}

calc_composition_prop_nt <- function(x, overlap = TRUE){
  
  n_counts <- x %>% letterFrequency(c("A", "C", "G", "T"))  # Nucleotide counts
  
  n_counts <- as_tibble(n_counts)
  n_counts$sums <- rowSums(n_counts)
  
  # normalise
  n_counts %<>% 
    dplyr::mutate(across(-sums),./sums) %>% 
    dplyr::select(-sums)
  
  df <- cbind(data.frame(
    nuc_id = x %>% names(),
    n_counts))  # Nucleotide counts
  return(df)
  
}

calc_composition_prop_dint <- function(x, overlap = TRUE){
  
  dn_counts_1 <- as_tibble(x %>%
                             DNAStringSet(start = 1) %>%
                             dinucleotideFrequency(step = 3) %>%
                             as.data.frame() %>%
                             rename_all(., ~ paste0(., "_p1")))
  
  dn_counts_1$sums <- rowSums(dn_counts_1)
  
  # normalise
  dn_counts_1 %<>% 
    dplyr::mutate(across(-sums),./sums) %>% 
    dplyr::select(-sums)
  
  dn_counts_2 <- as_tibble(x %>%
                             DNAStringSet(start = 2) %>%
                             dinucleotideFrequency(step = 3) %>%
                             as.data.frame() %>%
                             rename_all(., ~ paste0(., "_p2")))
  
  dn_counts_2$sums <- rowSums(dn_counts_2)
  
  # normalise
  dn_counts_2 %<>% 
    dplyr::mutate(across(-sums),./sums) %>% 
    dplyr::select(-sums)
  
  dn_counts_3 <- as_tibble(x %>%
                             DNAStringSet(start = 3) %>%
                             dinucleotideFrequency(step = 3) %>%
                             as.data.frame() %>%
                             rename_all(., ~ paste0(., "_p3")))
  
  dn_counts_3$sums <- rowSums(dn_counts_3)
  
  # normalise
  dn_counts_3 %<>% 
    dplyr::mutate(across(-sums),./sums) %>% 
    dplyr::select(-sums)
  
  df <- cbind(data.frame(
    nuc_id = x %>% names(),
    dn_counts_1, # Dinucleotide counts between positions 1-2 only
    dn_counts_2, # Dinucleotide counts between positions 2-3 only
    dn_counts_3)) # Dinucleotide counts between positions 3-1 only
  return(df)
  
}

calc_composition_prop_3mer <- function(x){
  
  kmers <- as_tibble(x %>% 
                       oligonucleotideFrequency(3, step=1))
  
  kmers$sums <- rowSums(kmers)
  
  # normalise
  kmers %<>% 
    dplyr::mutate(across(-sums),./sums) %>% 
    dplyr::select(-sums)
  
  df <- cbind(data.frame(
    nuc_id = x %>% names(),
    kmers))
  return(df)
  
}

calc_composition_prop_4mer <- function(x){
  
  kmers <- as_tibble(x %>% 
                       oligonucleotideFrequency(4, step=1))
  
  kmers$sums <- rowSums(kmers)
  
  # normalise
  kmers %<>% 
    dplyr::mutate(across(-sums),./sums) %>% 
    dplyr::select(-sums)
  
  df <- cbind(data.frame(
    nuc_id = x %>% names(),
    kmers))
  return(df)
  
}

calc_composition_prop_5mer <- function(x){
  
  kmers <- as_tibble(x %>% 
                       oligonucleotideFrequency(5, step=1))
  
  kmers$sums <- rowSums(kmers)
  
  # normalise
  kmers %<>% 
    dplyr::mutate(across(-sums),./sums) %>% 
    dplyr::select(-sums)
  
  df <- cbind(data.frame(
    nuc_id = x %>% names(),
    kmers))
  return(df)
  
}

calc_composition_prop_6mer <- function(x){
  
  kmers <- as_tibble(x %>% 
                       oligonucleotideFrequency(6, step=1))
  
  kmers$sums <- rowSums(kmers)
  
  # normalise
  kmers %<>% 
    dplyr::mutate(across(-sums),./sums) %>% 
    dplyr::select(-sums)
  
  df <- cbind(data.frame(
    nuc_id = x %>% names(),
    kmers))
  return(df)
  
}

### PROTEIN FEATS ###

calc_composition_prop_aac <- function(x){
  
  df <- cbind(data.frame(
    prot_id = x %>% names(),
    t(sapply(x, extractAAC)) # Calculate amino acid composition
  ))
  
  return(df)
}

calc_composition_prop_dp <- function(x){
  
  df <- cbind(data.frame(
    prot_id = x %>% names(),
    t(sapply(x, extractDC)) # Calculate dipeptide composition
  ))
  
  return(df)
}

calc_moranac_lag30 <- function(x){
  
  df <- cbind(data.frame(
    prot_id = x %>% names(),
    t(sapply(x, extractMoran,nlag=30))
  ))
  
  return(df)
}

calc_expanded_moranac_lag30 <- function(x){
  
  df <- cbind(data.frame(
    prot_id = x %>% names(),
    t(sapply(x, extractMoran,props=c("PRAM900102","PRAM900103","PRAM900104","CIDH920105","BHAR880101","CHAM820101","CHAM820102","CHOC760101",
                                     "BIGC670101","CHAM810101","JOND920102"),nlag=30)) # expanded set including frequency in secondary structural motifs and updated mutability index
  ))
  
  return(df)
}

calc_moranac_lag100 <- function(x){
  
  df <- cbind(data.frame(
    prot_id = x %>% names(),
    t(sapply(x, extractMoran,nlag=100))
  ))
  
  return(df)
}

calc_ctdc <- function(x){
  
  df <- cbind(data.frame(
    prot_id = x %>% names(),
    t(sapply(x, extractCTDC)) # Calculate CTD composition
  ))
  
  return(df)
}

calc_ctdt <- function(x){
  
  df <- cbind(data.frame(
    prot_id = x %>% names(),
    t(sapply(x, extractCTDT)) # Calculate CTD transition
  ))
  
  return(df)
}

calc_ctdd <- function(x){
  
  df <- cbind(data.frame(
    prot_id = x %>% names(),
    t(sapply(x, extractCTDD)) # Calculate CTD distribution
  ))
  
  return(df)
}

calc_ctriad <- function(x){
  
  df <- cbind(data.frame(
    prot_id = x %>% names(),
    t(sapply(x, extractCTriad)) # Calculate conjoint triads
  ))
  
  return(df)
}

### PLM EMBEDDING FEATS ###

read_embeddings <- function(x,n_pcomps){
  temp <- as_tibble(read_csv(x,col_names=F))
  
  # save the id then remove it from tibble
  id <- temp[1,ncol(temp)]
  temp <- temp[,-ncol(temp)]
  
  #drop the start and end tokens
  temp <- temp[-1,]
  temp <- temp[-nrow(temp),]
  
  temp %<>% dplyr::mutate_if(is.character,as.numeric) # convert all to numeric
  
  pc <- prcomp(temp,center =T,scale. =T)
  
  output <- tibble(principal_component = NA,pc_value = NA, ref_id = NA, seq_position = NA)
  
  for(i in 1:nrow(temp)){
    out <- tibble(principal_component = seq(1:n_pcomps),pc_value = pc$x[i,1:n_pcomps], ref_id = id$X961, seq_position = i)
    
    output <- rbind(output,out)
  }
  
  output <- output[-1,] # drop empty row
  
  return(output) # return the top n principal components
}

autocorrelate <- function(protein,pc,nlag){
  
  temp_df <- pc_embeddings %>% 
    dplyr::filter(ref_id == protein) %>% 
    dplyr::select(any_of(pc)) %>% 
    drop_na()
  
  temp_ac <- acf(temp_df,lag.max = nlag,plot=F)
  
  out <- tibble(lag = seq(1:nlag),ac_value = temp_ac$acf[1:nlag+1],principal_component = pc,ref_id = protein)   # drop the nlag 0 column, 
  
  return(out)  
}
