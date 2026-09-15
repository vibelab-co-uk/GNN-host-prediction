# define a function to apply function to every combination of 2 arguments
cmapply <- function(FUN, ..., MoreArgs = NULL){  
  # expand a grid of all argument combinations  
  l <- expand.grid(..., stringsAsFactors=FALSE)    
  # apply the function  
  .mapply(FUN=FUN, dots=unname(l), MoreArgs = MoreArgs)}

write_res <- function(results_dir){
  setwd(results_dir)
  
  res_filename <- list.files(pattern="^loo_out.*")
  
  # read in full results table
  result_preds <- read_csv(res_filename)
  
  # filter for only the prediction nodes
  result_preds %<>%
    filter(Train == F) # note that FALSE in the mask indicates a test node, therefore the prediction nodes are indicated by FALSE
  
  result_preds %<>%
    mutate(Correct = case_when(Observed == Predicted ~ T,
                               .default = F))
  
  result_preds$feat <- str_split(str_split(results_dir,pattern = "/")[[1]][3],pattern = "_")[[1]][2]
  
  write_csv(result_preds,paste0("loo_results_",result_preds$Hidden_layers[1],"_",result_preds$Model_structure[1],".csv"))
  
  setwd("../..")
}
