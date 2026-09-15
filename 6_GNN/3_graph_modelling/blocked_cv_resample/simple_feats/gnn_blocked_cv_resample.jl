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

include("./gnn_funcs.jl") # load loss and train functions

# get names of all feature sets with a regex
rx_file = r"count_feats.*"
rx_dir = r"results_.*"

in_feats = readdir("./feature_sets/", join=true)[occursin.(rx_file, readdir("./feature_sets/", join=true))]
out_dir = readdir("./feature_sets/", join=true)[occursin.(rx_dir, readdir("./feature_sets/", join=true))]

rx_blocks = r"reintro_.*"
resample_blocks = readdir("../../../1_generate_validation_sets/", join=true)[occursin.(rx_blocks, readdir("../../../1_generate_validation_sets/", join=true))]

# read in table of best tune values for each feature
best_tune = CSV.File("../../../../6_GNN/3_graph_modelling/kfold_cv_early_stopping/simple_feats/results_summary/best_tune_acc.csv") |> DataFrame

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
    infotime = 50      # report every `infotime` epochs
    nlayers = 2        # number of hidden layers
    patience = 60      # patience period before early stopping test applies

end

### read graph and metadata
# edgelist
coo_edgelist = readdlm("../../../../5_phylogeny/3_generate_graph/graph_data/patristic_paraymxo_graph.edgelist", Int64)

# GNN package doesn't allow 0 indexed edgelists
source = coo_edgelist[:, 1] .+ 1
target = coo_edgelist[:, 2] .+ 1

# edge weights
edge_weights = CSV.File("../../../../5_phylogeny/3_generate_graph/graph_data/patristic_paramyxo_graph_branch_weights.csv") |> DataFrame

# node metadata
node_names = CSV.File("../../../../5_phylogeny/3_generate_graph/graph_data/patristic_paramyxo_graph_node_names.csv") |> DataFrame
rename!(node_names, ["Node", "Accession"])

# host labels
host_labels = CSV.File("../../../../3_host_assignment/all_cluster_hosts.csv") |> DataFrame
rename!(host_labels, ["Accession", "Host_label"])

host_labels = dropmissing(leftjoin(node_names, host_labels, on=:Accession)) # use a join to ensure the host labels are in the same order as the node_names
sort!(host_labels, :Node)

# construct graph with edge weights
g = to_bidirected(GNNGraph(source, target, edge_weights.x |> f32))

# this gives a vector of host labels for each node
g.ndata.accessions = host_labels[:, 2]
g.ndata.labels = host_labels[:, 3]

# set up loop over blocked validation folds
for block in resample_blocks

    cv_folds = CSV.read(block, DataFrame, missingstring = "NA")

    function feature_comp_train(in_feats,out_dir,nhidden,nlayers)

        print("Training models with $nlayers hidden layers with dimension size $nhidden")

        node_feats = CSV.File(in_feats) |> DataFrame # read in node features
        node_feats = node_feats[:, Not(:Host_rank)] # remove host rank column from feature set

        node_feats = leftjoin(node_names, node_feats, on=:Accession) # use a join to ensure the node features are in the same order as the node_names
        sort!(node_feats, :Node)

        # convert the dataframe to an array with each column as feature vector for a node and assign node features
        g.ndata.features = transpose(Matrix{Float32}(node_feats[:, 3:end]))

        cv_out = DataFrame()# create empty DF to fill with ablated cv results

        # set up blocked cross validation loop across each subgraph
        for c in eachcol(cv_folds)
            # generate a boolean vector with true in the positions we want to use to train
            mask = BitVector(undef, nrow(node_feats)) # create an empty training mask
            mask .= true
            
            @. mask[c[!ismissing(c)]] = false # assign the test mask from the current CV fold
        
            # use the mask to assign train and test
            g.ndata.train_mask = mask
            g.ndata.test_mask = @. !mask

            ### train models
            abl_out = abl_train(g, out_dir, nhidden = nhidden, nlayers = nlayers)

            append!(cv_out,abl_out)
        end
        # add a column reporting which subgraph was blocked
        cv_out.cv_fold = repeat(names(cv_folds), inner = nrow(node_feats))
        
        resample_iter = filter.(isdigit, block) # extract the reintroduction count and iteration from the current filename

        cv_out.reintroductions .= resample_iter[2] # add these as ID columns and to the results file name
        cv_out.iteration .= resample_iter[3]

        CSV.write(joinpath(out_dir,"abl_output_reintro_$(resample_iter[2])_iter_$(resample_iter[3]).csv"),cv_out) # write the full ablation prediction results table
    end
    feature_comp_train.(in_feats, out_dir, best_tune.nhidden, best_tune.nlayers)
end
