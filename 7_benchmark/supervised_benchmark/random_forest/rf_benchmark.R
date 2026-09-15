library(tidyverse)
library(caret)
library(readxl)
library(writexl)
library(ranger)
library(magrittr)
library(doParallel)

### script to train a random forest model on multiple feature sets
`%nin%` <- Negate(`%in%`)

# set number of threads
numCoresAllowed <- 8

# open multithread cluster
cl <- makePSOCKcluster(numCoresAllowed)
registerDoParallel(cl)

# set seed for consistency
set.seed(4681)

# read in cross validation folds
kfolds <- read_csv("../../../6_GNN/1_generate_validation_sets/cv_folds_5.csv")

# define function to read in feature set and run model
feat_train <- function(feat_set){
  features <- read_csv(paste0("../feature_sets/",feat_set))
  
  host_ranks <- features %>% select(Accession,Host_rank)# create an ID/label vector to bind back to results
  
  feat_name <- str_sub(feat_set,start = 13, end = -5)
  
  num_feats <- ncol(features)-1 # number of features for tuning grid
  
  # # generate holdout validation set
  # inTraining <- createDataPartition(features$Host_rank, p = .8, list = FALSE)
  # train_set <- features[inTraining,]
  # test_set <- features[-inTraining,]
  # 
  # train_ids <- train_set$Accession # record order of accessions and drop the id column
  # train_set %<>% select(-Accession)
  # 
  # test_ids <- test_set$Accession # record order of accessions and drop the id column
  # test_set %<>% select(-Accession)
  
  # create cross validation folds
  cvFolds <- list()
  
  id_vec <- 1:nrow(host_ranks)
  
  id_vec[id_vec %nin% kfolds$Fold1]
  
  cvFolds$Fold1.Rep1 <- id_vec[id_vec %nin% kfolds$Fold1]
  cvFolds$Fold2.Rep1 <- id_vec[id_vec %nin% kfolds$Fold2]
  cvFolds$Fold3.Rep1 <- id_vec[id_vec %nin% kfolds$Fold3]
  cvFolds$Fold4.Rep1 <- id_vec[id_vec %nin% kfolds$Fold4]
  cvFolds$Fold5.Rep1 <- id_vec[id_vec %nin% kfolds$Fold5]
  
  features %<>%
    select(-Accession)
  
  # define train control with precision/recall summary function
  fitControl_cv <- trainControl(
    method="repeatedcv",
    classProbs = TRUE,
    verboseIter = TRUE,
    savePredictions = TRUE,
    returnResamp = "all",
    index = cvFolds
  )
  
  # set up tuning grid with mtry based on total number of features
  mtry_tune <- unique(c(round(num_feats/20),round(num_feats/10),round(num_feats/5),round(num_feats/2)))
  
  tgrid <- expand.grid(
    .mtry = mtry_tune,
    .splitrule = c("gini","extratrees"),
    .min.node.size = c(1,5,10)
  )
  
  # train model with cross validation on training set
  model_rf <- train(Host_rank ~ ., data=features,
                    method="ranger",
                    tuneGrid=tgrid,
                    trControl=fitControl_cv,
                    num.threads=1,
                    na.action=na.omit
  )
  
  # add column with the feature set used and write the training results
  train_res <- model_rf$results
  train_res$feature <- feat_name
  
  write.csv(train_res,paste0("results_",feat_name,"/rf_train_results.csv"))
  
  # # predict on test set
  # preds <- predict(model_rf,test_set,type="prob")
  # preds$Predicted <- predict(model_rf,test_set)
  # 
  # # add column of accession and observed host label
  # preds$Accession <- test_ids
  # 
  # preds <- left_join(preds,host_ranks)
  # 
  # write.csv(train_res,paste0("random_forest/results_",feat_name,"/rf_test_results.csv"))
}

dir_list <- list.files(path = "../feature_sets",pattern = "count_feats*")

lapply(dir_list,feat_train)

# close cluster
stopCluster(cl)
