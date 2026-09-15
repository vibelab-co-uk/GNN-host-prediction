# define a function to apply function to every combination of 2 arguments
cmapply <- function(FUN, ..., MoreArgs = NULL){  
  # expand a grid of all argument combinations  
  l <- expand.grid(..., stringsAsFactors=FALSE)    
  # apply the function  
  .mapply(FUN=FUN, dots=unname(l), MoreArgs = MoreArgs)}

combine_res <- function(dir){
  
  cv_list <- list()
  
  for(i in 1:max(cv_folds)){
    cv_list[[i]] <- read_csv(paste0("feature_sets/",dir,"/dumpout_output_cvfold",i,".csv"))
  }
  
  cv_res <- do.call(rbind,cv_list)
  
  write_csv(cv_res,paste0("feature_sets/",dir,"/dumpout_output_all.csv"))
  
  # extract the host label names
  host_labels <- str_sub(colnames(cv_res)[grepl("tax",colnames(cv_res))],start = 5)
  return(host_labels)
}

process_logs <- function(file,dir){
  # read in training log and extract data
  train_log <- read_delim(paste0("./feature_sets/",dir,"/",file),delim="\\n",col_names = F)
  
  test_log <- str_split_i(train_log$X1,pattern = "Test",i = 2)
  train_log <- str_split_i(train_log$X1,pattern = "Test",i = 1)
  
  epoch = unlist(str_extract_all(train_log,pattern = "(?<=Epoch\\:\\s)\\d+"))
  train_loss = unlist(str_extract_all(train_log,pattern = "(?<=loss\\s\\=\\s)\\d+\\.\\d+"))
  train_accuracy = unlist(str_extract_all(train_log,pattern = "(?<=acc\\s\\=\\s)\\d+\\.\\d+"))
  
  test_loss = unlist(str_extract_all(test_log,pattern = "(?<=loss\\s\\=\\s)\\d+\\.\\d+"))
  test_accuracy = unlist(str_extract_all(test_log,pattern = "(?<=acc\\s\\=\\s)\\d+\\.\\d+"))
  
  # make a tibble
  full_log <- tibble("epoch" = epoch,"training_loss" = train_loss, "training_acc" = train_accuracy,"test_loss" = test_loss, "test_acc" = test_accuracy)
  
  full_log$feat <- str_sub(dir,start=9)
  
  full_log$cv_fold <- str_extract_all(file,"\\d+")[[1]][length(str_extract_all(file,"\\d+")[[1]])]
  return(full_log)
}

# define a function to read in model training logs
read_acc_loss <- function(files,dir){

  full_log <- do.call(rbind,lapply(files,FUN = process_logs,dir = dir))
  
  return(full_log)
}

process_tsne <- function(file,dir){
  # read in tsne table
  tsne <- read_csv(paste0("./feature_sets/",dir,"/",file),col_names=F)
  
  tsne$cv_fold <- str_extract_all(file,"\\d+")[[1]][length(str_extract_all(file,"\\d+")[[1]])]
  
  tsne <- cbind(tsne,tip_hosts)
  
  tsne$epoch <- str_extract_all(file,"\\d+")[[1]][length(str_extract_all(file,"\\d+")[[1]])-1]
  return(tsne)
}

# define a function to read in t-sne data
read_tsne <- function(files,dir){
  
  feat <- str_sub(dir,start=9)
  
  all_tsne <- do.call(rbind,lapply(files,FUN = process_tsne,dir = dir))
  
  all_tsne$feat <- feat
  
  return(all_tsne)
}

tsne_plots <- function(x){
  
  ggplot(x,aes(x=tSNE1,y=tSNE2,col=Observed)) +
    geom_point(alpha=0.5)+
    scale_color_manual(values=pal)+
    theme_bw(base_size = 12)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    facet_wrap(~epoch)
  
  ggsave(paste0("./results_summary/tsne_plots/tsne_plot_",x$cv_fold[1],"_",x$feat[1],".png"),height = 1920, width = 2880,units = "px")
}


# define function to read in raw embedding tables and calculate the distance between raw embeddings and patristic distances
calc_distances <- function(file,dir,hosts){

  raw_embed <- read_csv(paste0("./feature_sets/",dir,"/",file),col_names=F)

  raw_embed <- as.matrix(t(raw_embed))
  
  rownames(raw_embed) <- tip_hosts$ref_accessions
  
  # return euclidean distance between all output embeddings
  embed_dist <- as.matrix(dist(raw_embed,method = "euclidean"))
  
  # normalise the distance matrices
  embed_dist <- (embed_dist-min(embed_dist))/(max(embed_dist)-min(embed_dist))
  
  # calculate residuals & norm
  sum_resid <- sum(embed_dist-pat_distances)
  
  sum_square_resid <- sum((embed_dist-pat_distances) ^ 2)
  
  L2_norm <- sqrt(sum((embed_dist-pat_distances) ^ 2))
  
  # # calculate the correlation using the mantel test
  # mantel_cor <- mantel(pat_distances,embed_dist,method = "spearman",permutations = 999)
  # 
  # distance_tbl <- tibble("sum_residuals" = sum(embed_dist-pat_distances),"sum_square_residuals" = sum((embed_dist-pat_distances) ^ 2),
  #                        "L2_norm" = sqrt(sum((embed_dist-pat_distances) ^ 2)),"mantel_correlation" = mantel_cor$statistic,
  #                        "mantel_significance" = mantel_cor$signif, epoch = str_extract_all(file,"\\d+")[[1]][length(str_extract_all(file,"\\d+")[[1]])-1],
  #                        cv_fold = str_extract_all(file,"\\d+")[[1]][length(str_extract_all(file,"\\d+")[[1]])],feature_set=dir)
  
  # short version not including mantel statistics which greatly increase runtime
  distance_tbl <- tibble("sum_residuals" = sum(embed_dist-pat_distances),"sum_square_residuals" = sum((embed_dist-pat_distances) ^ 2),
                         "L2_norm" = sqrt(sum((embed_dist-pat_distances) ^ 2)),epoch = str_extract_all(file,"\\d+")[[1]][length(str_extract_all(file,"\\d+")[[1]])-1],
                         cv_fold = str_extract_all(file,"\\d+")[[1]][length(str_extract_all(file,"\\d+")[[1]])],feature_set=str_sub(dir,start=9))
  
  
  return(distance_tbl)
}

# define function to calculate distance metrics across all conditions and return dataframes for plotting
pat_embedding_similarity <- function(files,dir){
  
  # calculate distances
  all_distance <- do.call(rbind,lapply(files,FUN = calc_distances,dir = dir,hosts = tip_hosts))
  
  return(all_distance)
}

plot_max_min_results <- function(dir){
  
  feat <- str_sub(dir,start=9)
  
  all_res <- read_csv(paste0("./feature_sets/",dir,"/dumpout_output_all.csv"))
  
  feat_dist <- all_embedding_dist %>% 
    filter(feature_set == feat)
  
  # find the epochs where patristic distance and embeddings are least similar
  min_dist <- feat_dist %>% 
    filter(sum_square_residuals == min(sum_square_residuals,na.rm = T))
  
  max_dist <- feat_dist %>% 
    filter(sum_square_residuals == max(sum_square_residuals))
  
  min_res <- all_res %>% 
    filter(Training_epochs == min_dist$epoch)
  
  max_res <- all_res %>% 
    filter(Training_epochs == max_dist$epoch)
  
  min_preds <- min_res %>%
    group_by(Accession) %>% 
    summarise(across(starts_with("tax"),~mean(.x)))  
  
  min_preds %<>% 
    group_by(Accession) %>% 
    mutate(Prediction_score = max(across(starts_with("tax")))) %>% 
    # subsequent mutate assigns the name of the tax column with the max score (thereby giving the prediction)
    mutate(Predicted = str_sub(colnames(min_preds)[grepl("^tax",colnames(min_preds))][which.max(across(starts_with("tax")))],start = 5))
  
  # bind back in some metadata
  min_preds <- left_join(min_preds,min_res %>% select(Accession,Observed) %>% distinct())
  
  min_preds %<>%
    mutate(outcome = case_when(Predicted == Observed ~ "Correct",
                               .default = "Incorrect"))
  
  ### and max
  max_preds <- max_res %>%
    group_by(Accession) %>% 
    summarise(across(starts_with("tax"),~mean(.x)))  
  
  max_preds %<>% 
    group_by(Accession) %>% 
    mutate(Prediction_score = max(across(starts_with("tax")))) %>% 
    # subsequent mutate assigns the name of the tax column with the max score (thereby giving the prediction)
    mutate(Predicted = str_sub(colnames(max_preds)[grepl("^tax",colnames(max_preds))][which.max(across(starts_with("tax")))],start = 5))
  
  # bind back in some metadata
  max_preds <- left_join(max_preds,max_res %>% select(Accession,Observed) %>% distinct())
  
  max_preds %<>%
    mutate(outcome = case_when(Predicted == Observed ~ "Correct",
                               .default = "Incorrect"))
  
  
  
  minshape <- shapes[factor(min_preds$outcome,levels = c("Correct","Incorrect"))]
  maxshape <- shapes[factor(max_preds$outcome,levels = c("Correct","Incorrect"))]
  
  png(filename=paste0("./results_summary/max_min_residuals_dist_plots/min_graph_results_",feat,".png"),width=1200,height=2400,units="px")
  plot(g, vertex.size=3,vertex.color=mycol,vertex.shape = minshape,vertex.label=NA,edge.width=3,layout = layout)
  legend('topright',legend=levels(factor((node_id$Host_rank))),col = palette,pch = 16,cex=3)
  dev.off()
  
  png(filename=paste0("./results_summary/max_min_residuals_dist_plots/max_graph_results_",feat,".png"),width=1200,height=2400,units="px")
  plot(g, vertex.size=3,vertex.color=mycol,vertex.shape = maxshape,vertex.label=NA,edge.width=3,layout = layout)
  legend('topright',legend=levels(factor((node_id$Host_rank))),col = palette,pch = 16,cex=3)
  dev.off()
}

plot_correct_label_transitions <- function(dir){
  
  feat <- str_sub(dir,start=9)
  
  all_res <- read_csv(paste0("./feature_sets/",dir,"/dumpout_output_all.csv"))
  
  labels <- all_res %>% select(Accession,Observed) %>% distinct()
  
  # calculate results across cv folds
  all_res %<>%
    group_by(Accession,Training_epochs) %>% 
    summarise(across(starts_with("tax"),~mean(.x)))  
  
  all_res %<>% 
    group_by(Accession,Training_epochs) %>% 
    mutate(Prediction_score = max(across(starts_with("tax")))) %>% 
    # subsequent mutate assigns the name of the tax column with the max score (thereby giving the prediction)
    mutate(Predicted = str_sub(colnames(all_res)[grepl("^tax",colnames(all_res))][which.max(across(starts_with("tax")))],start = 5)) %>% 
    select(Accession,Training_epochs,Prediction_score,Predicted)
  
  # bind back in ground truth labels
  all_res <- left_join(all_res,labels)
  
  all_res %<>% 
    mutate(Outcome = case_when(Predicted == Observed ~ "Correct",
                               .default = "Incorrect"))
  
  all_res$Outcome <- factor(all_res$Outcome,levels = c("Incorrect","Correct"))
  
  all_res$Observed <- factor(all_res$Observed,levels = c("Aves","Rodentia","Chiroptera","Artiodactyla","Eulipotyphla","Carnivora",
                                                         "Primates","Reptilia","Fish","Perissodactyla","Pholidota","Scandentia")) # rerank according to frequency
  
  ggplot(data = all_res,aes(x = Training_epochs,y = Outcome,group = Accession))+
    geom_jitter(width = 0,alpha = 0.2)+
    geom_line(linetype = "dashed",color = "firebrick")+
    theme_bw() +
    xlab("Training epochs") +
    facet_wrap(~Observed) +
    theme(strip.background = element_rect(fill="white"))
  
  ggsave(paste0("./results_summary/correctness_scatters/correctness_barcode_scatter_",feat,".png"))
  
  # in order to optimise the next step, filter for only those predictions which change during the final 25% of model training
  transition_res <- all_res %>% 
    group_by(Accession) %>% 
    filter(Training_epochs >= quantile(Training_epochs)[4]) %>% 
    mutate(transition = length(unique(Outcome))) %>% 
    filter(transition != 1)
  
  switches <- do.call(rbind,lapply(unique(transition_res$Accession),label_transitions, transitions = transition_res))
  
  switches %<>%
    select(-transition) %>% 
    filter(transition_switch %in% c("change_in_error","correct_change"))
  
  switches$feat <- feat
  message(paste0("Feature ",feat," processed succesfully"))
  return(switches)
}

label_transitions <- function(trans_accession,transitions){
  
  target_res <- transitions %>% 
    filter(Accession == trans_accession)
  
  target_res$transition_switch <- NA
  
  for(i in 2:nrow(target_res)){
    if(target_res$Outcome[i-1] == target_res$Outcome[i]){
      target_res$transition_switch[i] <- "no_change"
    } else if (target_res$Outcome[i-1] == "Correct" & target_res$Outcome[i] == "Incorrect") {
      target_res$transition_switch[i] <- "change_in_error"
    } else {
      target_res$transition_switch[i] <- "correct_change"
    }
  }
  
  return(target_res)
}
