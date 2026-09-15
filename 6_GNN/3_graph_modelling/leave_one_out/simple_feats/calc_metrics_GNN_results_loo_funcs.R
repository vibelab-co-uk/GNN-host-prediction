# define a function to apply function to every combination of 2 arguments
cmapply <- function(FUN, ..., MoreArgs = NULL){  
  # expand a grid of all argument combinations  
  l <- expand.grid(..., stringsAsFactors=FALSE)    
  # apply the function  
  .mapply(FUN=FUN, dots=unname(l), MoreArgs = MoreArgs)}

write_res <- function(results_dir){
  setwd(results_dir)
  
  res_filename <- list.files(pattern="^ensemble_out.*")
  
  # read in full results table
  result_preds <- read_csv(res_filename)
  
  # filter for only the prediction nodes
  result_preds %<>%
    filter(Prediction_node == F) # note that FALSE in the mask indicates a test node, therefore the prediction nodes are indicated by FALSE
  
  # make final prediction across ensemble folds
  ensemble_preds <- result_preds %>%
    group_by(Accession) %>% 
    summarise(across(starts_with("tax"),~mean(.x))) %>% 
    group_by(Accession) %>% 
    mutate(Prediction_score = max(across(starts_with("tax")))) %>% 
    # subsequent mutate assigns the name of the tax column with the max score (thereby giving the prediction)
    mutate(Predicted = str_sub(colnames(result_preds)[grepl("^tax",colnames(result_preds))][which.max(across(starts_with("tax")))],start = 5))
  
  # bind back in some metadata
  ensemble_preds <- left_join(ensemble_preds,result_preds %>% select(Accession,Observed,Hidden_layers,Model_structure) %>% distinct())

  ensemble_preds %<>%
    mutate(Correct = case_when(Observed == Predicted ~ T,
                               .default = F))
  
  ensemble_preds$feat <- str_split(str_split(results_dir,pattern = "/")[[1]][3],pattern = "_")[[1]][2]
  
  write_csv(ensemble_preds,paste0("ensemble_results_",result_preds$Hidden_layers[1],"_",result_preds$Model_structure[1],".csv"))

  setwd("../..")
}
