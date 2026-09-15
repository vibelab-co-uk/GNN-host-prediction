# EVAL LOSS FUNCTION
function eval_loss_accuracy(X, y, mask, model)
    ŷ = model(X)
    l = logitcrossentropy(ŷ[:, mask], y[:, mask])
    acc = mean(onecold(ŷ[:, mask]) .== onecold(y[:, mask]))
    return (loss = round(l, digits = 4), acc = round(acc * 100, digits = 2))
end

function nn_dumpout_train(features, cv_fold, truth, out_dir, mask, accessions; kws...)
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

    X = features
    y = onehotbatch(truth |> cpu, unique(truth)) |> device
    ytrain = y[:, mask]

    classes = unique(truth) # make a unique list of classes to name outputs

    nin, nout = size(X, 1), length(unique(truth))

    if args.nlayers == 1
        model = Chain(Dense(nin => args.nhidden, relu),
                     Dense(Int64(args.nhidden/2) => nout)) |> device
     elseif args.nlayers == 2
        model = Chain(Dense(nin => args.nhidden, relu),
                     Dense(args.nhidden => Int64(args.nhidden/2),relu),
                     Dense(Int64(args.nhidden/2) => nout)) |> device
    elseif args.nlayers == 3
        model = Chain(Dense(nin => args.nhidden, relu),
                     Dense(args.nhidden => Int64(args.nhidden/2),relu),
                     Dense(Int64(args.nhidden/2) => Int64(args.nhidden/4),relu),
                     Dense(Int64(args.nhidden/4) => nout)) |> device
     elseif args.nlayers == 4
        model = Chain(Dense(nin => args.nhidden, relu),
                     Dense(args.nhidden => Int64(args.nhidden/2),relu),
                     Dense(Int64(args.nhidden/2) => Int64(args.nhidden/4),relu),
                     Dense(Int64(args.nhidden/4) => Int64(args.nhidden/8),relu),
                     Dense(Int64(args.nhidden/8) => nout)) |> device
    end

    opt = Flux.setup(Adam(args.η), model)

    # LOGGING FUNCTION
    function report(epoch)
        train = eval_loss_accuracy(X, y, mask, model)
        test = eval_loss_accuracy(X, y, .!mask, model)
        open(joinpath(out_dir,"training_log_$(args.nhidden)_$(args.nlayers)convlayers_weighted.txt"),"a") do io
            println(io,"Epoch: $epoch   Train: $(train)   Test: $(test)")
        end
    end

    # EARLY STOPPING FUNCTION
    function terminate(stop_vec)
        early_stop = false
        if length(stop_vec) >= 10
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
    z = tsne(model(X) |> transpose,2) 
    CSV.write(joinpath(out_dir,"raw_embeddings_$(args.nhidden)_$(args.nlayers)_epoch0_cvfold$cv_fold.csv"),DataFrame(model(X),:auto),append = true)
    CSV.write(joinpath(out_dir,"tsne_$(args.nhidden)_$(args.nlayers)_epoch0_cvfold$cv_fold.csv"),DataFrame(z,:auto),append = true)

    ## assess trained model
    final_out = DataFrame(Observed = truth, Prediction_score = onecold(model(X)), Predicted = Vector{String}(undef,length(truth)),
                        tax1 = softmax(model(X))[1, :], tax2 = softmax(model(X))[2, :], tax3 = softmax(model(X))[3, :],
                        tax4 = softmax(model(X))[4, :], tax5 = softmax(model(X))[5, :], tax6 = softmax(model(X))[6, :],
                        tax7 = softmax(model(X))[7, :], tax8 = softmax(model(X))[8, :], tax9 = softmax(model(X))[9, :],
                        tax10 = softmax(model(X))[10, :], tax11 = softmax(model(X))[11, :], tax12 = softmax(model(X))[12, :],
                        Train = mask, Hidden_layers = args.nhidden, Model_structure = string(args.nlayers)*"denselayers", Accession = accessions, 
                        Training_epochs = length(stop_vec))

    tax_names = "tax_" .* classes

    rename!(final_out, ["tax$i" => tax for (i, tax) in enumerate(tax_names)])

    # we have R code to calculate micro/macro CM stats. Build tables with probabilities and preds in Julia then write the output and do downstream analysis in R
    final_out.Prediction_score = select(final_out,r"tax" => ByRow(max) => :Prediction_score).Prediction_score

    # calculate predicted class for every row by finding the column that matches the max value
    for i in names(select(final_out,r"tax")), l in 1:nrow(final_out)
        final_out.Prediction_score[l] == final_out[!, i][l] ? final_out.Predicted[l] = i[5:end] : missing   
    end
    # finally add a column reporting which cv fold results came from
    final_out.cv_fold .= cv_fold

    CSV.write(joinpath(out_dir,"dumpout_output_cvfold$cv_fold.csv"),final_out) # write the full ablation prediction results table

    for epoch in 1:(args.epochs)
        grad = Flux.gradient(model) do model
            ŷ = model(X)
            logitcrossentropy(ŷ[:, mask], ytrain)
        end

        Flux.update!(opt, model, grad[1])

        # append current model loss for the whole graph to the stopping vector
        append!(stop_vec, eval_loss_accuracy(X, y, trues(length(truth)), model)[1,])

        # on report epochs write tsne
        if epoch % args.infotime == 0 
            z = tsne(model(X) |> transpose,2) 
            CSV.write(joinpath(out_dir,"tsne_$(args.nhidden)_$(args.nlayers)_epoch$(epoch)_cvfold$cv_fold.csv"),DataFrame(z,:auto),append = true)
        end

        # dump out raw embeddings and training log with loss and accuracy every epoch
        CSV.write(joinpath(out_dir,"raw_embeddings_$(args.nhidden)_$(args.nlayers)_epoch$(epoch)_cvfold$cv_fold.csv"),DataFrame(model(X),:auto),append = true)
        
        report(epoch)

        ## assess trained model
        final_out = DataFrame(Observed = truth, Prediction_score = onecold(model(X)), Predicted = Vector{String}(undef,length(truth)),
                            tax1 = softmax(model(X))[1, :], tax2 = softmax(model(X))[2, :], tax3 = softmax(model(X))[3, :],
                            tax4 = softmax(model(X))[4, :], tax5 = softmax(model(X))[5, :], tax6 = softmax(model(X))[6, :],
                            tax7 = softmax(model(X))[7, :], tax8 = softmax(model(X))[8, :], tax9 = softmax(model(X))[9, :],
                            tax10 = softmax(model(X))[10, :], tax11 = softmax(model(X))[11, :], tax12 = softmax(model(X))[12, :],
                            Train = mask, Hidden_layers = args.nhidden, Model_structure = string(args.nlayers)*"denselayers", Accession = accessions, 
                            Training_epochs = length(stop_vec))

        tax_names = "tax_" .* classes

        rename!(final_out, ["tax$i" => tax for (i, tax) in enumerate(tax_names)])

        # we have R code to calculate micro/macro CM stats. Build tables with probabilities and preds in Julia then write the output and do downstream analysis in R
        final_out.Prediction_score = select(final_out,r"tax" => ByRow(max) => :Prediction_score).Prediction_score

        # calculate predicted class for every row by finding the column that matches the max value
        for i in names(select(final_out,r"tax")), l in 1:nrow(final_out)
            final_out.Prediction_score[l] == final_out[!, i][l] ? final_out.Predicted[l] = i[5:end] : missing   
        end
        
        # finally add a column reporting which cv fold results came from
        final_out.cv_fold .= cv_fold

        CSV.write(joinpath(out_dir,"dumpout_output_cvfold$cv_fold.csv"),final_out,append=true) # write the full ablation prediction results table

        # test for early stopping
        if length(stop_vec) >= args.patience && terminate(stop_vec) == true # if patience period has expired and early stopping function returns true, stop training
            @goto escape_label
        end
    end

    # if early stopping criterion is met the training loop exits early
    @label escape_label

end