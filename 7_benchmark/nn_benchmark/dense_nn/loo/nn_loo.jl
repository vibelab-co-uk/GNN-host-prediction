begin
    using CUDA
    using MLDatasets, MLUtils
    using LinearAlgebra, Random, Statistics, StatsBase, GLM
    using DelimitedFiles, DataFrames, CSV
    using Flux
    using Flux: onehotbatch, onecold, logitcrossentropy
end

include("./nn_funcs.jl") # load loss and train functions

## get names of all feature sets with a regex
rx_file = r"count_feats.*"
rx_dir = r"results_.*"

in_feats = readdir("./feature_sets/", join=true)[occursin.(rx_file, readdir("./feature_sets/", join=true))]
out_dir = readdir("./feature_sets/", join=true)[occursin.(rx_dir, readdir("./feature_sets/", join=true))]

# read in table of best tune values for each feature
best_tune = CSV.File("../kfold_cv/results_summary/best_tune_acc.csv") |> DataFrame

# sort to match input feats/out directories (alphabetical)
sort!(best_tune, [order(:feat)])

# update dataframe for 5 random pseudofeatures
best_tune = DataFrame(nhidden = [best_tune.nhidden;[best_tune.nhidden[nrow(best_tune)],best_tune.nhidden[nrow(best_tune)],best_tune.nhidden[nrow(best_tune)],best_tune.nhidden[nrow(best_tune)]]],
                      nlayers = [best_tune.nlayers;[best_tune.nlayers[nrow(best_tune)],best_tune.nlayers[nrow(best_tune)],best_tune.nlayers[nrow(best_tune)],best_tune.nlayers[nrow(best_tune)]]],
                      feat = [best_tune.feat;[best_tune.feat[nrow(best_tune)],best_tune.feat[nrow(best_tune)],best_tune.feat[nrow(best_tune)],best_tune.feat[nrow(best_tune)]]])

# define struct with arguments for the `train` function 
Base.@kwdef mutable struct Args
    η = 1.0f-3             # learning rate
    epochs = 500          # number of epochs
    seed = 22             # set seed > 0 for reproducibility
    nhidden = 128        # dimension size of hidden layers
    nlayers = 2          # number of hidden layers
    usecuda = true      # if true use cuda (if available)
    infotime = 50      # report every `infotime` epochs
    patience = 60      # patience period before early stopping test applies
    pdropout = 0.2     # proportion dropout
end
##
function feature_comp_train(in_feats,out_dir,nhidden,nlayers)

    feats = CSV.File(in_feats) |> DataFrame # read in node features
    host_labels = feats.Host_rank
    accessions = feats.Accession
    feats = feats[:, Not(:Host_rank,:Accession)] # remove host rank and accession column from feature set

    # convert the dataframe to an array with each column as feature vector
    feats = transpose(Matrix{Float32}(feats))

    loo_out = DataFrame()# create empty DF to fill with ablated cv results

    # set up k-fold cross validation loop
    for i in 1:length(host_labels)
        # generate a boolean vector with true in the positions we want to use to train
        mask = BitVector(undef, length(host_labels)) # create an empty training mask
        mask .= true
        
        mask[i] = false # assign the prediction mask from the current prediction node

        ### train models
        nn_out = nn_train(feats, host_labels, out_dir, mask, accessions, nhidden = nhidden, nlayers = nlayers)

        append!(loo_out,nn_out)
    end
    
    CSV.write(joinpath(out_dir,"nn_loo_output.csv"),loo_out) # write the full ablation prediction results table
end
feature_comp_train.(in_feats, out_dir, best_tune.nhidden, best_tune.nlayers)
