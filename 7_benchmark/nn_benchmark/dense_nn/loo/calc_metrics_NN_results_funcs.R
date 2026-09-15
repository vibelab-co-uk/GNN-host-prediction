# define a function to apply function to every combination of 2 arguments
cmapply <- function(FUN, ..., MoreArgs = NULL){  
  # expand a grid of all argument combinations  
  l <- expand.grid(..., stringsAsFactors=FALSE)    
  # apply the function  
  .mapply(FUN=FUN, dots=unname(l), MoreArgs = MoreArgs)}

# assign functions for performance metrics --------------------------------
# accuracy
calculate.accuracy <- function(x) {
  x$Predicted <- as.character(x$Predicted)
  x$Observed <- as.character(x$Observed)
  
  acc_vec <- length(which(x$Predicted == x$Observed)) / length(x$Observed)
  return(tibble("Accuracy"=acc_vec))
}

# calculate confusion matrix stats from cm
get.conf.stats <- function(cm) {
  out <- vector("list", length(cm))
  for (i in seq_along(cm)) {
    x <- cm[[i]]
    tp <- x$table[x$positive, x$positive]
    fp <- sum(x$table[x$positive, colnames(x$table) != x$positive])
    fn <- sum(x$table[colnames(x$table) != x$positive, x$positive])
    # TNs are not well-defined for one-vs-all approach
    elem <- c(tp = tp, fp = fp, fn = fn)
    out[[i]] <- elem
  }
  df <- do.call(rbind, out)
  rownames(df) <- unlist(lapply(cm, function(x) x$positive))
  return(as.data.frame(df))
}

# calculate micro f1
get.micro.f1 <- function(cm) {
  cm.summary <- get.conf.stats(cm)
  tp <- sum(cm.summary$tp)
  fn <- sum(cm.summary$fn)
  fp <- sum(cm.summary$fp)
  pr <- tp / (tp + fp)
  re <- tp / (tp + fn)
  f1 <- 2 * ((pr * re) / (pr + re))
  return(f1)
}

# calculate macro f1
get.macro.f1 <- function(cm) {
  c <- cm[[1]]$byClass # a single matrix is sufficient
  re <- sum(c[, "Recall"], na.rm=T) / nrow(c)
  pr <- sum(c[, "Precision"], na.rm=T) / nrow(c)
  f1 <- 2 * ((re * pr) / (re + pr))
  return(f1)
}

micro_macro_f1 <- function(protein){
  cm <- vector("list", length(levels(protein$Observed)))
  for (i in seq_along(cm)) {
    positive.class <- levels(protein$Observed)[i]
    # in the i-th iteration, use the i-th class as the positive class
    cm[[i]] <- confusionMatrix(protein$Predicted, protein$Observed, 
                               positive = positive.class)
    
  }
  micro.f1 <- get.micro.f1(cm)
  macro.f1 <- get.macro.f1(cm)
  
  return(tibble("micro_f1" = micro.f1,"macro_f1" = macro.f1,"protein" = unique(protein$Protein)))
}

factorise <- function(protein){
  protein$Predicted <- factor(protein$Predicted,levels=c("Aves","Rodentia","Chiroptera","Artiodactyla","Primates",
                                                         "Eulipotyphla","Carnivora","Perissodactyla","Scandentia",
                                                         "Pholidota","Fish","Reptilia"))
  
  protein$Observed <- factor(protein$Observed,levels=c("Aves","Rodentia","Chiroptera","Artiodactyla","Primates",
                                                       "Eulipotyphla","Carnivora","Perissodactyla","Scandentia",
                                                       "Pholidota","Fish","Reptilia"))
  return(protein)
}

OvA_ROC <- function(tbl,target){
  
  temp <- tbl %>% 
    mutate(Observed = case_when(Observed == target ~ "Positive",
                                .default = "Negative"))
  
  # set up mleval table
  target_prob <- temp %>% 
    dplyr::select(contains(as.character(target)))  
  
  target_prob <- as.vector(target_prob %>% unname())
  
  roc_tbl <- tibble("Positive"=target_prob[[1]],"Negative"=1-target_prob[[1]],"obs"=temp$Observed,"Group"=target)
  
  return(roc_tbl)
}

subgraph_analysis <- function(abl_tbl,subgraph){
  # unless we are handling the whole graph (subgraph 0), filter for target subgraph
  if(subgraph != 0){
    abl_tbl %<>%
      dplyr::filter(Subgraph == subgraph)
  }
  
  subgraph_tax <- unique(c(abl_tbl$Observed,abl_tbl$Predicted))
  
  # drop any groups for which there are no observations in this subgraph
  `%nin%` <- Negate(`%in%`)
  
  to_drop <- levels(subgraph_tax)[levels(subgraph_tax) %nin% subgraph_tax]
  
  sub_abl_tbl <- abl_tbl
  
  # reassign factor using only levels present in this subgraph
  sub_abl_tbl$Observed <- factor(sub_abl_tbl$Observed, levels = subgraph_tax)
  sub_abl_tbl$Predicted <- factor(sub_abl_tbl$Predicted, levels = subgraph_tax)
  
  # for cases where we have only a single class in either group, we cannot derive f1
  if(length(unique(sub_abl_tbl$Observed)) == 1 | length(unique(sub_abl_tbl$Predicted)) == 1){
    
    agg_results_order <- cbind(tibble(micro_f1 = NA, macro_f1 = NA, tax_level = "Order"),calculate.accuracy((sub_abl_tbl)))
    
  } else {
    
    agg_results_order <- cbind(micro_macro_f1(sub_abl_tbl),calculate.accuracy(sub_abl_tbl))
    
    agg_results_order$tax_level <- "Order"
  }
  
  # repeat this analysis at the level of taxonomic class
  result_preds_class <- sub_abl_tbl %>%
    mutate(Predicted = case_when(Predicted == "Aves" ~ factor("Aves"),
                                 Predicted == "Fish" ~ factor("Fish"),
                                 Predicted == "Reptilia" ~ factor("Reptilia"),
                                 .default = factor("Mammal")), 
           Observed = case_when(Observed == "Aves" ~ factor("Aves"),
                                Observed == "Fish" ~ factor("Fish"),
                                Observed == "Reptilia" ~ factor("Reptilia"),
                                .default = factor("Mammal")))
  
  # for cases where we have only a single class, we cannot derive f1
  if(length(unique(result_preds_class$Observed)) == 1 | length(unique(result_preds_class$Predicted)) == 1){
    
    agg_results_class <- cbind(tibble(micro_f1 = NA, macro_f1 = NA, tax_level = "Class"),calculate.accuracy((result_preds_class)))
    
  } else {
    
    agg_results_class <- cbind(micro_macro_f1(result_preds_class),calculate.accuracy(result_preds_class))
    agg_results_class$tax_level <- "Class"
  }
  
  agg_results <- rbind(agg_results_order,agg_results_class)
  
  agg_results$Ablation <- abl_tbl$drop_feat[1]
  agg_results$Subgraph <- subgraph
  
  rm(agg_results_order,agg_results_class)
  
  full_roc_tbl <- do.call(rbind,lapply(all_tax,OvA_ROC,tbl=abl_tbl))
  
  full_roc_tbl$obs <- factor(full_roc_tbl$obs, levels = c("Positive","Negative"))
  
  subgraph_tax <- unique(abl_tbl$Observed)
  
  # drop any groups for which there are no observations in this subgraph
  `%nin%` <- Negate(`%in%`)
  
  to_drop <- levels(subgraph_tax)[levels(subgraph_tax) %nin% subgraph_tax]
  
  full_roc_tbl %<>%
    filter(Group %nin% to_drop)
  
  # reassign factor using only levels present in this subgraph
  full_roc_tbl$Group <- factor(full_roc_tbl$Group, levels = subgraph_tax)
  
  # for cases where we have only a single observed class, we manually assign metric scores
  if(length(unique(full_roc_tbl$obs)) == 1 & length(unique(full_roc_tbl$Group) == 1)){
    if(unique(full_roc_tbl$obs) == "Positive"){
      full_out <- tibble(Metric = c("SENS", "SPEC", "MCC", "Informedness", "PREC", "NPV", "FPR", "F1", "TP", "FP",          
                                    "TN", "FN", "AUC-ROC", "AUC-PR", "AUC-PRG"),
                         Score = c(NA,NA,NA,NA,1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA), CI = rep(NA,times=15), 
                         Taxa = subgraph_tax, Subgraph = subgraph)
    } else {
      full_out <- tibble(Metric = c("SENS", "SPEC", "MCC", "Informedness", "PREC", "NPV", "FPR", "F1", "TP", "FP",          
                                    "TN", "FN", "AUC-ROC", "AUC-PR", "AUC-PRG"),
                         Score = c(NA,NA,NA,NA,0,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA), CI = rep(NA,times=15), 
                         Taxa = subgraph_tax, Subgraph = subgraph)
    }
    
    message(paste0("Subgraph ",subgraph,"; some metrics could not be calculated due to small sample size"))
  } else {
    
    full_eval <- evalm(data.frame(full_roc_tbl),plots='roc',rlinethick=0.8,fsize=16,bins=12,positive="Positive")
    
    full_out <- tibble(do.call(rbind, full_eval$stdres))
    
    full_out$metric <- rownames(do.call(rbind, full_eval$stdres))
    
    full_out <- separate_wider_delim(full_out,metric,delim=".",names=c("Taxa","Metric"))
    
    full_out$Subgraph <- subgraph
    
    full_out %<>%
      dplyr::select("Metric","Score","CI","Taxa","Subgraph")
    
    message(paste0("Subgraph ",subgraph,"; processed succesfully"))
  }
  
  ### calculate the no info rates ###
  # find most frequent class
  class_counts <- sub_abl_tbl %>% 
    dplyr::count(Observed)
  
  no_info_tbl <- sub_abl_tbl
  no_info_tbl$Predicted <- class_counts$Observed[which.max(class_counts$n)] # assign the most frequent class as predicted class in all cases
  
  # calculate metrics for the no info tbl
  # for cases where we have only a single class in either group, we cannot derive f1
  if(length(unique(no_info_tbl$Observed)) == 1 | length(unique(no_info_tbl$Predicted)) == 1){
    
    agg_no_info_order <- cbind(tibble(micro_f1 = NA, macro_f1 = NA, tax_level = "Order"),calculate.accuracy((no_info_tbl)))
    
  } else {
    
    agg_no_info_order <- cbind(micro_macro_f1(no_info_tbl),calculate.accuracy(no_info_tbl))
    
    agg_no_info_order$tax_level <- "Order"
  }
  
  # repeat this analysis at the level of taxonomic class
  # find most frequent class
  class_counts <- result_preds_class %>% 
    dplyr::count(Observed)
  
  sub_abl_tbl_class <- sub_abl_tbl %>%
    mutate(Predicted = case_when(Predicted == "Aves" ~ factor("Aves"),
                                 Predicted == "Fish" ~ factor("Fish"),
                                 Predicted == "Reptilia" ~ factor("Reptilia"),
                                 .default = factor("Mammal")), 
           Observed = case_when(Observed == "Aves" ~ factor("Aves"),
                                Observed == "Fish" ~ factor("Fish"),
                                Observed == "Reptilia" ~ factor("Reptilia"),
                                .default = factor("Mammal")))
  
  no_info_class <- sub_abl_tbl_class
  no_info_class$Predicted <- class_counts$Observed[which.max(class_counts$n)] # assign the most frequent class as predicted class in all cases
  
  # for cases where we have only a single class, we cannot derive f1
  if(length(unique(no_info_class$Observed)) == 1 | length(unique(no_info_class$Predicted)) == 1){
    
    agg_no_info_class <- cbind(tibble(micro_f1 = NA, macro_f1 = NA, tax_level = "Class"),calculate.accuracy((result_preds_class)))
    
  } else {
    
    agg_no_info_class <- cbind(micro_macro_f1(no_info_class),calculate.accuracy(no_info_class))
    agg_no_info_class$tax_level <- "Class"
  }
  
  agg_no_info <- rbind(agg_no_info_order,agg_no_info_class)
  
  agg_no_info$Subgraph <- subgraph
  
  return(list(agg_results,full_out,agg_no_info))
}

write_res <- function(results_dir){
  setwd(results_dir)
  
  master_acc_tbl <- tibble("Accuracy"=NA,"Protein"=NA,"tune_params"=NA)
  master_f1_tbl <- tibble("micro_f1"=NA,"macro_f1"=NA,"protein"=NA,"tune_params"=NA)
  
  # read in full results table
  res_filename <- list.files(pattern="^nn_out*")
  
  result_preds <- read_csv(res_filename)
  
  result_preds <- left_join(result_preds,subgraph_membership)
  
  result_preds <- factorise(result_preds)
  
  # lazy superassignment used to avoid scoping issues
  all_tax <<- unique(result_preds$Observed)
  
  # apply analysis function across all ablated features and subgraphs
  abl_out <- cmapply(FUN = subgraph_analysis,subgraph = c(0,unique(subgraph_membership$Subgraph)),MoreArgs = list(abl_tbl = result_preds))
  
  agg_results_full <- tibble()
  
  for(i in 1:length(abl_out)){
    agg_results_full <- rbind(agg_results_full,abl_out[[i]][[1]])
  }
  
  ROC_results_full <- tibble()
  
  for(i in 1:length(abl_out)){
    ROC_results_full <- rbind(ROC_results_full,abl_out[[i]][[2]])
  }
  
  noinfo_results_full <- tibble()
  
  for(i in 1:length(abl_out)){
    noinfo_results_full <- rbind(noinfo_results_full,abl_out[[i]][[3]])
  }
  
  write_csv(agg_results_full,paste0("aggregate_results.csv"))
  write_csv(ROC_results_full,paste0("ROC_PR_results.csv"))
  write_csv(noinfo_results_full,paste0("no_info_results.csv"))
  
  setwd("../..")
}

cat_agg <- function(results_dir){
  setwd(results_dir)
  
  agg_results <- read_csv(paste0("aggregate_results.csv"))
  
  agg_results$feat <- str_split(str_split(results_dir,pattern = "/")[[1]][3],pattern = "_")[[1]][2]

  setwd("../..")
  return(agg_results)
}

cat_roc <- function(results_dir){
  setwd(results_dir)
  
  full_out <- read_csv(paste0("ROC_PR_results.csv"))
  
  full_out$feat <- str_split(str_split(results_dir,pattern = "/")[[1]][3],pattern = "_")[[1]][2]

  setwd("../..")
  return(full_out)
}

cat_noinfo <- function(results_dir){
  setwd(results_dir)
  
  full_out <- read_csv(paste0("no_info_results.csv"))
  
  full_out$feat <- str_split(str_split(results_dir,pattern = "/")[[1]][3],pattern = "_")[[1]][2]

  setwd("../..")
  return(full_out)
}
