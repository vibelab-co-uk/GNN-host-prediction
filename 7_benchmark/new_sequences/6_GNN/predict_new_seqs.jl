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
    using CUDA
    using BSON
end

# get names of all feature sets with a regex
rx_file = r"count_feats.*"
rx_dir = r"results_.*"

in_feats = readdir("./feature_sets/", join=true)[occursin.(rx_file, readdir("./feature_sets/", join=true))]
out_dir = readdir(".", join=true)[occursin.(rx_dir, readdir(".", join=true))]

### read graph with new sequences and metadata
# edgelist
coo_edgelist = readdlm("../5_phylogeny/3_generate_graph/graph_data/patristic_paraymxo_graph.edgelist", Int64)

# GNN package doesn't allow 0 indexed edgelists
source = coo_edgelist[:, 1] .+ 1
target = coo_edgelist[:, 2] .+ 1

# edge weights
edge_weights = CSV.File("../5_phylogeny/3_generate_graph/graph_data/patristic_paramyxo_graph_branch_weights.csv") |> DataFrame

# node metadata
node_names = CSV.File("../5_phylogeny/3_generate_graph/graph_data/patristic_paramyxo_graph_node_names.csv") |> DataFrame
rename!(node_names, ["Node", "Accession"])

# host labels
host_labels = CSV.File("../3_host_assignment/combined_cluster_hosts.csv") |> DataFrame
rename!(host_labels, ["Accession", "Host_label"])

host_labels = dropmissing(leftjoin(node_names, host_labels, on=:Accession)) # use a join to ensure the host labels are in the same order as the node_names
sort!(host_labels, :Node)

# construct graph with edge weights
g = to_bidirected(GNNGraph(source, target, edge_weights.x |> f32))

# this gives a vector of host labels for each node
g.ndata.accessions = host_labels[:, 2]
g.ndata.labels = host_labels[:, 3]

# load in relevant model and dataset and predict on new sequences
function predict_new(in_feats, out_dir)
    
    # load model
    feature_id = chop(out_dir,head = 10,tail = 0)

    rx_bson = r"trained_model_.*"

    path = "../../../6_GNN/3_graph_modelling/full_data_train/simple_feats/feature_sets/results_"*feature_id

    bson = readdir(path, join=true)[occursin.(rx_bson, readdir(path, join=true))]

    BSON.@load bson[1] model

    # load features 
    node_feats = CSV.File(in_feats) |> DataFrame # read in node features
    node_feats = node_feats[:, Not(:Host_rank)] # remove host rank column from feature set

    node_feats = leftjoin(node_names, node_feats, on=:Accession) # use a join to ensure the node features are in the same order as the node_names
    sort!(node_feats, :Node)

    # convert the dataframe to an array with each column as feature vector for a node and assign node features
    g.ndata.features = transpose(Matrix{Float32}(node_feats[:, 3:end]))

    # get model hyperparameters from filename
    model_params = split(bson[1],"/")
    model_params = model_params[length(model_params)]

    nhidden = match(r"(?<=nhidden)\d+",model_params).match
    nlayers = match(r"(?<=nlayers)\d+",model_params).match
    pdropout = match(r"(?<=dropout)\d+\.?\d*",model_params).match
    epochs = match(r"(?<=epochs)\d+",model_params).match

    # make predictions
    X = g.features

    # assign classes, using same order as in model training
    classes = ["Primates", "Rodentia", "Carnivora", "Scandentia", "Chiroptera", "Reptilia", "Fish", "Artiodactyla", "Aves", "Perissodactyla", "Eulipotyphla", "Pholidota"]    

    final_out = DataFrame(Observed = g.labels, Prediction_score = onecold(model(g, X)), Predicted = Vector{String}(undef,length(g.labels)),
                        tax1 = softmax(model(g, X))[1, :], tax2 = softmax(model(g, X))[2, :], tax3 = softmax(model(g, X))[3, :],
                        tax4 = softmax(model(g, X))[4, :], tax5 = softmax(model(g, X))[5, :], tax6 = softmax(model(g, X))[6, :],
                        tax7 = softmax(model(g, X))[7, :], tax8 = softmax(model(g, X))[8, :], tax9 = softmax(model(g, X))[9, :],
                        tax10 = softmax(model(g, X))[10, :], tax11 = softmax(model(g, X))[11, :], tax12 = softmax(model(g, X))[12, :],
                        Hidden_layers = nhidden, Model_structure = nlayers, Accession = g.ndata.accessions, 
                        Training_epochs = epochs, Dropout = pdropout)

    tax_names = "tax_" .* classes

    rename!(final_out, ["tax$i" => tax for (i, tax) in enumerate(tax_names)])

    # we have R code to calculate micro/macro CM stats. Build tables with probabilities and preds in Julia then write the output and do downstream analysis in R
    final_out.Prediction_score = select(final_out,r"tax" => ByRow(max) => :Prediction_score).Prediction_score

    # calculate predicted class for every row by finding the column that matches the max value
    for i in names(select(final_out,r"tax")), l in 1:nrow(final_out)
        final_out.Prediction_score[l] == final_out[!, i][l] ? final_out.Predicted[l] = i[5:end] : missing   
    end
    CSV.write(joinpath(out_dir,"new_seqs_output.csv"),final_out) # write the full ablation prediction results table

end

predict_new.(in_feats, out_dir)
