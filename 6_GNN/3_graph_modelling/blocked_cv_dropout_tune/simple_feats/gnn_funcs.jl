function eval_loss_accuracy(X, y, mask, model, g)
    ŷ = model(g, X)
    l = logitcrossentropy(ŷ[:, mask], y[:, mask])
    acc = mean(onecold(ŷ[:, mask]) .== onecold(y[:, mask]))
    return (loss = round(l, digits = 4), acc = round(acc * 100, digits = 2))
end

function abl_train(g, out_dir; kws...)
    args = Args(; kws...)

    args.seed > 0 && Random.seed!(args.seed)

    if args.usecuda && CUDA.functional()
        device = gpu
        args.seed > 0 && CUDA.seed!(args.seed)
        @info "Training on GPU"
    else
        device = cpu
        @info "Training on CPU"
    end

    # set up a vector to fill with loss values to use for early stopping test
    stop_vec = Vector{Float64}()
    
    X = g.features
    y = onehotbatch(g.labels |> cpu, unique(g.ndata.labels)) |> device
    ytrain = y[:, g.ndata.train_mask]

    classes = unique(g.ndata.labels) # make a unique list of classes to name outputs

    nin, nhidden, nout = size(X, 1), args.nhidden, length(unique(g.ndata.labels))

    # define the GNN model based on the nlayers argument
    if args.nlayers == 1
        model = GNNChain(GCNConv(nin => nhidden, relu, use_edge_weight = true),
                        Dropout(args.pdropout),
                        Dense(nhidden, nout)) |> device
    elseif args.nlayers == 2
        model = GNNChain(GCNConv(nin => nhidden, relu, use_edge_weight = true),
                        GCNConv(nhidden => nhidden, relu, use_edge_weight = true),
                        Dropout(args.pdropout),
                        Dense(nhidden, nout)) |> device
    elseif args.nlayers == 3
        model = GNNChain(GCNConv(nin => nhidden, relu, use_edge_weight = true),
                        GCNConv(nhidden => nhidden, relu, use_edge_weight = true),
                        GCNConv(nhidden => nhidden, relu, use_edge_weight = true),
                        Dropout(args.pdropout),
                        Dense(nhidden, nout)) |> device
    elseif args.nlayers == 4
        model = GNNChain(GCNConv(nin => nhidden, relu, use_edge_weight = true),
                        GCNConv(nhidden => nhidden, relu, use_edge_weight = true),
                        GCNConv(nhidden => nhidden, relu, use_edge_weight = true),
                        GCNConv(nhidden => nhidden, relu, use_edge_weight = true),
                        Dropout(args.pdropout),
                        Dense(nhidden, nout)) |> device
    else 
        println("ERROR: This number of layers is not supported")
    end

    opt = Flux.setup(Adam(args.η), model)

    display(g)

    ## LOGGING FUNCTION
    function report(epoch)
        train = eval_loss_accuracy(X, y, g.train_mask, model, g)
        test = eval_loss_accuracy(X, y, g.test_mask, model, g)
        open(joinpath(out_dir,"training_log_$(args.pdropout).txt"),"a") do io
            println(io,"Epoch: $epoch   Train: $(train)   Test: $(test)")
        end
    end

    ## EARLY STOPPING FUNCTION
    function terminate(stop_vec)
        early_stop = false
        if args.nhidden >= 512 && args.nlayers >= 3 && length(stop_vec) >= 30 # for larger models, we have more stringent stopping criteria
            Y = stop_vec[length(stop_vec)-29:length(stop_vec)]

            stop_df = DataFrame(X = Float64.(1:length(Y)), Y = Y)
            stop_lm = lm(@formula(Y ~ X),stop_df) # calculate rate of change of loss over last 30 epochs

            loss_rate = GLM.coef(stop_lm)[2] - stderror(stop_lm)[2]

            # early stop if loss rate is >= 0 over last 30 epochs and 10 consecutive epochs with no decrease in model loss
            if loss_rate >= 0 && 
                stop_vec[length(stop_vec)] >= stop_vec[length(stop_vec)-1] && 
                stop_vec[length(stop_vec)-1] >= stop_vec[length(stop_vec)-2] && 
                stop_vec[length(stop_vec)-2] >= stop_vec[length(stop_vec)-3] && 
                stop_vec[length(stop_vec)-3] >= stop_vec[length(stop_vec)-4] && 
                stop_vec[length(stop_vec)-4] >= stop_vec[length(stop_vec)-5] &&
                stop_vec[length(stop_vec)-5] >= stop_vec[length(stop_vec)-6] &&
                stop_vec[length(stop_vec)-6] >= stop_vec[length(stop_vec)-7] &&
                stop_vec[length(stop_vec)-7] >= stop_vec[length(stop_vec)-8] &&
                stop_vec[length(stop_vec)-8] >= stop_vec[length(stop_vec)-9] &&
                stop_vec[length(stop_vec)-9] >= stop_vec[length(stop_vec)-10]

                early_stop = true
            end
        elseif length(stop_vec) >= 10
            Y = stop_vec[length(stop_vec)-9:length(stop_vec)]

            stop_df = DataFrame(X = Float64.(1:length(Y)), Y = Y)
            stop_lm = lm(@formula(Y ~ X),stop_df) # calculate rate of change of loss over last 10 epochs
            
            loss_rate = GLM.coef(stop_lm)[2] - stderror(stop_lm)[2]
            
            # early stop if loss rate is >= 0 over last 10 epochs
            if loss_rate >= 0 
                early_stop = true
            end
        end
        early_stop
    end

    ## TRAINING
    report(0)
    for epoch in 1:(args.epochs)
        grad = Flux.gradient(model) do model
            ŷ = model(g, X)
            logitcrossentropy(ŷ[:, g.train_mask], ytrain)
        end

        Flux.update!(opt, model, grad[1])

        # append current model loss for the whole graph to the stopping vector
        append!(stop_vec, eval_loss_accuracy(X, y, trues(212), model, g)[1,])
        
        # report
        epoch % args.infotime == 0 && report(epoch)

        # test for early stopping
        if length(stop_vec) >= args.patience && terminate(stop_vec) == true # if patience period has expired and early stopping function returns true, stop training
            @goto escape_label
        end
    end

    # if early stopping criterion is met the training loop exits early
    @label escape_label

    ## assess trained model
    final_out = DataFrame(Observed = g.labels, Prediction_score = onecold(model(g, X)), Predicted = Vector{String}(undef,length(g.labels)),
                        tax1 = softmax(model(g, X))[1, :], tax2 = softmax(model(g, X))[2, :], tax3 = softmax(model(g, X))[3, :],
                        tax4 = softmax(model(g, X))[4, :], tax5 = softmax(model(g, X))[5, :], tax6 = softmax(model(g, X))[6, :],
                        tax7 = softmax(model(g, X))[7, :], tax8 = softmax(model(g, X))[8, :], tax9 = softmax(model(g, X))[9, :],
                        tax10 = softmax(model(g, X))[10, :], tax11 = softmax(model(g, X))[11, :], tax12 = softmax(model(g, X))[12, :],
                        Train = g.train_mask, Hidden_layers = nhidden, Model_structure = string(args.nlayers)*"convlayers_weighted", 
                        Accession = g.ndata.accessions, Training_epochs = length(stop_vec))

    tax_names = "tax_" .* classes

    rename!(final_out, ["tax$i" => tax for (i, tax) in enumerate(tax_names)])

    # we have R code to calculate micro/macro CM stats. Build tables with probabilities and preds in Julia then write the output and do downstream analysis in R
    final_out.Prediction_score = select(final_out,r"tax" => ByRow(max) => :Prediction_score).Prediction_score

    # calculate predicted class for every row by finding the column that matches the max value
    for i in names(select(final_out,r"tax")), l in 1:nrow(final_out)
        final_out.Prediction_score[l] == final_out[!, i][l] ? final_out.Predicted[l] = i[5:end] : missing   
    end
    final_out
end