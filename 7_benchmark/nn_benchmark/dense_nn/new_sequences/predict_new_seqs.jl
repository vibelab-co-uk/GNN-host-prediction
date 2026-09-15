begin
    using CUDA
    using MLDatasets, MLUtils
    using LinearAlgebra, Random, Statistics, StatsBase, GLM
    using DelimitedFiles, DataFrames, CSV
    using Flux
    using Flux: onehotbatch, onecold, logitcrossentropy
    using BSON
end

# get names of all feature sets with a regex
rx_file = r"count_feats.*"
rx_dir = r"results_.*"

in_feats = readdir("./feature_sets/", join=true)[occursin.(rx_file, readdir("./feature_sets/", join=true))]
out_dir = readdir(".", join=true)[occursin.(rx_dir, readdir(".", join=true))]

# load in relevant model and dataset and predict on new sequences
function predict_new(in_feats, out_dir)
    
    # load model
    feature_id = chop(out_dir,head = 10,tail = 0)

    rx_bson = r"trained_model_.*"

    path = "../full_data_train/feature_sets/results_"*feature_id

    bson = readdir(path, join=true)[occursin.(rx_bson, readdir(path, join=true))]

    BSON.@load bson[1] model

    # load features 
    feats = CSV.File(in_feats) |> DataFrame # read in node features

    host_labels = feats.Host_rank
    accessions = feats.Accession

    feats = feats[:, Not(:Host_rank,:Accession)] # remove host rank and accession column from feature set

    # convert the dataframe to an array with each column as feature vector for a node and assign node features
    feats = transpose(Matrix{Float32}(feats))

    # get model hyperparameters from filename
    model_params = split(bson[1],"/")
    model_params = model_params[length(model_params)]

    nhidden = match(r"(?<=nhidden)\d+",model_params).match
    nlayers = match(r"(?<=nlayers)\d+",model_params).match
    pdropout = match(r"(?<=dropout)\d+\.?\d*",model_params).match
    epochs = match(r"(?<=epochs)\d+",model_params).match

    # make predictions
    X = feats

    # assign classes, using same order as in model training
    classes = ["Primates", "Rodentia", "Carnivora", "Scandentia", "Chiroptera", "Reptilia", "Fish", "Artiodactyla", "Aves", "Perissodactyla", "Eulipotyphla", "Pholidota"]    

    ## assess trained model
    final_out = DataFrame(Observed = host_labels, Prediction_score = onecold(model(X)), Predicted = Vector{String}(undef,length(host_labels)),
                        tax1 = softmax(model(X))[1, :], tax2 = softmax(model(X))[2, :], tax3 = softmax(model(X))[3, :],
                        tax4 = softmax(model(X))[4, :], tax5 = softmax(model(X))[5, :], tax6 = softmax(model(X))[6, :],
                        tax7 = softmax(model(X))[7, :], tax8 = softmax(model(X))[8, :], tax9 = softmax(model(X))[9, :],
                        tax10 = softmax(model(X))[10, :], tax11 = softmax(model(X))[11, :], tax12 = softmax(model(X))[12, :],
                        Hidden_layers = nhidden, Model_structure = nlayers, Accession = accessions, 
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
