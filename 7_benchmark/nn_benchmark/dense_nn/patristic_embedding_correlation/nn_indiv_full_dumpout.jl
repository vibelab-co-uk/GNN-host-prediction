begin
    using Graphs
    using GraphIO.EdgeList
    using DataFrames, CSV
    using GLM, StatsBase
    using Flux
    using Flux: onecold, onehotbatch, logitcrossentropy
    using LinearAlgebra, Random, Statistics, StatsBase
    using GraphNeuralNetworks
    using DelimitedFiles
    using MLDatasets, MLUtils
    using TSne, Distances
    using CUDA
end

include("./nn_funcs_full_dumpout.jl") # load loss and train functions

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
                      feat = [best_tune.feat;[best_tune.feat[nrow(best_tune)],best_tune.feat[nrow(best_tune)],best_tune.feat[nrow(best_tune)],best_tune.feat[nrow(best_tune)]]])

# define struct with arguments for the `train` function 
Base.@kwdef mutable struct Args
    η = 1.0f-3             # learning rate
    epochs = 500          # number of epochs
    seed = 22             # set seed > 0 for reproducibility
    usecuda = true      # if true use cuda (if available)
    nhidden = 128       # dimension of hidden features
    infotime = 25      # report every `infotime` epochs
    nlayers = 2        # number of hidden layers
    patience = 60      # patience period before early stopping test applies
end

# read in graph metadata
# this will be used to arrange features to ensure CV folds are the same as for the graph models
# node metadata
node_names = CSV.File("../../../../5_phylogeny/3_generate_graph/graph_data/patristic_paramyxo_graph_node_names.csv") |> DataFrame
rename!(node_names, ["Node", "Accession"])

# graph node host labels
node_hosts = CSV.File("../../../../3_host_assignment/all_cluster_hosts.csv") |> DataFrame
rename!(node_hosts, ["Accession", "Host_label"])

node_hosts = dropmissing(leftjoin(node_names, node_hosts, on=:Accession)) # use a join to ensure the host labels are in the same order as the node_names
sort!(node_hosts, :Node)

# read in cross validation folds
cv_folds = CSV.read("../../../../6_GNN/1_generate_validation_sets/cv_folds_5.csv", DataFrame, missingstring = "NA")

function feature_comp_train(in_feats,out_dir,nhidden,nlayers)

    feats = CSV.File(in_feats) |> DataFrame # read in node features

    feats = leftjoin(node_hosts, feats, on=:Accession) # use a join to ensure same ordering as for the graph models
    sort!(feats, :Node)

    host_labels = feats.Host_label
    accessions = feats.Accession

    feats = feats[:, Not(:Host_rank,:Host_label,:Accession,:Node)] # remove host rank and accession column from feature set

    # convert the dataframe to an array with each column as feature vector
    feats = transpose(Matrix{Float32}(feats))

    # set up k-fold cross validation loop
    for c in 1:ncol(cv_folds)
        # generate a boolean vector with true in the positions we want to use to train
        mask = BitVector(undef, length(host_labels)) # create an empty training mask
        mask .= true
        
        current_fold = cv_folds[:,c]

        @. mask[current_fold[!ismissing(current_fold)]] = false # assign the test mask from the current CV fold

        ### train models
        nn_out = nn_dumpout_train(feats, c, host_labels, out_dir, mask, accessions, nhidden = nhidden, nlayers = nlayers)
    end
end

feature_comp_train.(in_feats, out_dir, best_tune.nhidden, best_tune.nlayers)
