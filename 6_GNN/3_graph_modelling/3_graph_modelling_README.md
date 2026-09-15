Process for training and performance assessment on GNN models\
Each subdirectory trains models with different feature sets or validation strategies\
All use the same three basic scripts\
1: run the Julia script to train and output model results e.g. `gnn_indiv_kfold_cv.jl` or `gnn_blocked_cv.jl`\
2: run the R script which calculates machine learning performance metrics and saves results e.g. `calc_metrics_GNN_results.R` or `calc_metrics_GNN_results_blocked_test_only`\
3: run the R script which generates figures and summary tables e.g. `analyse_GNN_results.R`}

The `patristic_embedding_correlation` directory only has a single analysis script `embedding_analysis.R` which carries out all analysis and generates figures