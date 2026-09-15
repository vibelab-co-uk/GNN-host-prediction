function eval_loss_accuracy(X, y, model)
    ŷ = model(X)
    l = logitcrossentropy(ŷ, y)
    acc = mean(onecold(ŷ) .== onecold(y))
    return (loss = round(l, digits = 4), acc = round(acc * 100, digits = 2))
end

function nn_train(features, truth, out_dir, accessions; kws...)
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

    classes = unique(truth) # make a unique list of classes to name outputs

    nin, nout = size(X, 1), length(unique(truth))

    if args.nlayers == 1
        model = Chain(Dense(nin => args.nhidden, relu),
                     Dropout(args.pdropout),
                     Dense(Int64(args.nhidden/2) => nout)) |> device
     elseif args.nlayers == 2
        model = Chain(Dense(nin => args.nhidden, relu),
                     Dense(args.nhidden => Int64(args.nhidden/2),relu),
                     Dropout(args.pdropout),
                     Dense(Int64(args.nhidden/2) => nout)) |> device
    elseif args.nlayers == 3
        model = Chain(Dense(nin => args.nhidden, relu),
                     Dense(args.nhidden => Int64(args.nhidden/2),relu),
                     Dense(Int64(args.nhidden/2) => Int64(args.nhidden/4),relu),
                     Dropout(args.pdropout),
                     Dense(Int64(args.nhidden/4) => nout)) |> device
     elseif args.nlayers == 4
        model = Chain(Dense(nin => args.nhidden, relu),
                     Dense(args.nhidden => Int64(args.nhidden/2),relu),
                     Dense(Int64(args.nhidden/2) => Int64(args.nhidden/4),relu),
                     Dense(Int64(args.nhidden/4) => Int64(args.nhidden/8),relu),
                     Dropout(args.pdropout),
                     Dense(Int64(args.nhidden/8) => nout)) |> device
    end

    opt = Flux.setup(Adam(args.η), model)

    ## LOGGING FUNCTION
    function report(epoch)
        train = eval_loss_accuracy(X, y, model)
        open(joinpath(out_dir,"training_log_$(args.nhidden)_$(args.nlayers)_$(args.pdropout)_convlayers_weighted.txt"),"a") do io
            println(io,"Epoch: $epoch   Train: $(train)")
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
            ŷ = model(X)
            logitcrossentropy(ŷ, y)
        end

        Flux.update!(opt, model, grad[1])

        # append current model loss for the whole graph to the stopping vector
        append!(stop_vec, eval_loss_accuracy(X, y, model)[1,])
        
        # report
        epoch % args.infotime == 0 && report(epoch)

        # test for early stopping
        if length(stop_vec) >= args.patience && terminate(stop_vec) == true # if patience period has expired and early stopping function returns true, stop training
            @goto escape_label
        end
    end
    
    # if early stopping criterion is met the training loop exits early
    @label escape_label

    @save joinpath(out_dir,"trained_model_nhidden$(args.nhidden)_nlayers$(args.nlayers)_dropout$(args.pdropout)_epochs$(length(stop_vec)).bson") model
end