begin
    using CUDA
    using MLDatasets, MLUtils
    using LinearAlgebra, Random, Statistics, StatsBase, GLM
    using DelimitedFiles, DataFrames, CSV
    using Flux
    using Flux: onehotbatch, onecold, logitcrossentropy
    using BSON: @save
end

include("./nn_funcs_save_model.jl") # load loss and train functions

# get names of all feature sets with a regex
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
                      feat = [best_tune.feat;[best_tune.feat[nrow(best_tune)],best_tune.feat[nrow(best_tune)],best_tune.feat[nrow(best_tune)],best_tune.feat[nrow(best_tune)]]],
                      epochs = [best_tune.epochs;[best_tune.epochs[nrow(best_tune)],best_tune.epochs[nrow(best_tune)],best_tune.epochs[nrow(best_tune)],best_tune.epochs[nrow(best_tune)]]])

# define struct with arguments for the `train` function 
Base.@kwdef mutable struct Args
    η = 1.0f-3             # learning rate
    epochs = 500          # number of epochs
    seed = 22             # set seed > 0 for reproducibility
    usecuda = true      # if true use cuda (if available)
    nhidden = 128       # dimension of hidden features
    infotime = 50      # report every `infotime` epochs
    nlayers = 2        # number of hidden layers
    patience = 60      # patience period before early stopping test applies
    pdropout = 0.2     # proportion dropout 
end


function feature_comp_train(in_feats,out_dir,nhidden,nlayers,epochs)

    feats = CSV.File(in_feats) |> DataFrame # read in node features

    host_labels = feats.Host_rank
    accessions = feats.Accession

    feats = feats[:, Not(:Host_rank,:Accession)] # remove host rank and accession column from feature set

    # convert the dataframe to an array with each column as feature vector
    feats = transpose(Matrix{Float32}(feats))

    ### train models
    nn_out = nn_train(feats, host_labels, out_dir, accessions, nhidden = nhidden, nlayers = nlayers, epochs = epochs)

end

feature_comp_train.(in_feats, out_dir, best_tune.nhidden, best_tune.nlayers, best_tune.epochs)
