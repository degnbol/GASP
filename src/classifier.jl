#!/usr/bin/env julia
using ArgParse
using DataFrames, CSV
using Statistics
using Chain: @chain
using Statistics, Random
using Flux
using Flux: logitbinarycrossentropy
ROOT = readchomp(`git root`)
include("$ROOT/src/utilities/glob.jl")
include("$ROOT/src/utilities/AUC.jl")
include("$ROOT/src/utilities/dataframe.jl")


"""
Custom metric to evaluate predictive performance.
Has biotechnological relevance, e.g. for predicting top P candidates, and testing them in lab.
Metric is on a scale from 0 to 1, where 1 is perfect prediction on top P datapoints, 
and 0 means no true positives among top P datapoints.
Approach:
1.	Set P equal to number of positive datapoints.
2.	Sort datapoints by prediction score.
3.	Assign weights 1…P to top P datapoints from lowest prediction score to highest.
4.	Normalize weights so they sum to 1.
5.	The performance metric is the sum of weights for positive datapoints.
"""
function topP_metric(pred::Vector{<:AbstractFloat}, class::Union{BitVector,Vector{Bool}})::Float64
    P = sum(class)
    class = class[sortperm(pred)][end-P+1:end]
    weights = (1:P) ./ sum(1:P)
    weights[class] |> sum
end
# unit test
@assert(topP_metric(
    [0.12, 0.3, 0.31, 0.4, 0.5, 0.6, 0.7, 0.72, 0.8, 0.9], 
    [false, false, false, true, true, false, true, false, true, false]) == 0.4,
    "Unit test failed")

parser = ArgParseSettings(description=
"""Train (-o/--train) or evaluate (-i/--eval) a random forest classifier.""")
@add_arg_table parser begin
    "infiles"
    arg_type = String # default is Any
    nargs = '+'
    help = "Infile with a datapoint on each row for training or prediction,
    the former requiring the presence of a class column (-c/--class).
    Extra files supplied should be feature annotation files,
    which will be (left)joined onto the first (left)file using any shared columns given by -j/--join."
    "--delim", "-d"
    arg_type = Char
    default = '\t'
    "--join", "-j"
    arg_type = String
    nargs = '+'
    help = "Delimiter separating columns in table files."
    "--eval", "-i"
    arg_type = String
    help = "Evaluation: read previously trained model from path."
    "--train", "-o"
    arg_type = String
    help = "Training: train model and write it to path."
    "--class", "-c"
    arg_type = String
    nargs = '+'
    help = "Required for training. Column name containing true class.
    Multiple can be given for getting cor and AUC on testset, 
    where the first column provided found in trainset is used as the training class."
    "--model", "-m"
    arg_type = String
    nargs = '+'
    default = ["n_trees=100", "partial_sampling=1."]
    help = "Provide args given to ScikitLearn RandomForestClassifier.
    https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html"
    "--test", "-t"
    help = "For quick evaluation, i.e. train and prediction in one call.
    If a path is given it is assumed to be a table file similar to the first infile.
    If a string is given found as a column name in any of the infiles,
    then the column should contain 0 and 1, where 1 indicates testset,
    or (not implemented yet) more than two unique values which each indicate sets for cross-validation.
    If an integer/fraction is given then this many/fraction of datapoints will be randomly chosen for testing."
    "--substitute", "-s"
    arg_type = String
    nargs = 2
    help = "Encode a column with a substitution matrix.
    First arg is name of column to encode,
    second arg is path to substitution matrix,
    either as absolute path, relative to pwd or to data/NCBI/.
    Should contain numeric columns for substitution and a single string column for naming the new features.
    Default file extension is .tsv
    so simply \"blosum62Amb\" or \"matchAmb\" will be enough."
end

# if run as script
if abspath(PROGRAM_FILE) == @__FILE__
    args = ARGS
else
    # if interactive, do an example
    ROOT = readchomp(`git root`)
    pTrainset = glob("$ROOT/results/*-cmpPerform/RF/trainset.tsv") |> only
    pTestset = glob("$ROOT/results/*-cmpPerform/RF/testset.tsv") |> only
    pChems = glob("$ROOT/results/*-chemicalFeatures/acceptors2_features.tsv") |> only
    sArgs = "$pTrainset $pChems -j cid enzyme -t $pTestset -c reaction rate -s seq blosum62Amb -o RF.jld2"
    args = split(sArgs, ' ')
end
args = parse_args(args, parser, as_symbols=true) |> NamedTuple

# check args
@assert (length(args.infiles) < 2) || (args.join !== nothing) "Multiple infiles requires -j/--join"
@assert sum([args.train, args.eval] .!== nothing) == 1 "Provide exactly one of -i/--train or -o/--eval"
@assert (args.train === nothing) || (args.class !== nothing) "-c/--class is required for training"

readdf(source::AbstractString) = CSV.read(source, DataFrame; delim=args.delim)

# edge case where a column contains only missing entries
dropmissing_cols(df::DataFrame) = df[!, eltype.(eachcol(df)) .!= Missing]

function dropmissing_verbose!(df::DataFrame, cols=:)
    cols isa Colon || intersect!(cols, names(df))
    colHasMiss = [n for n in names(df[!, cols]) if any(ismissing, df[!, n])]
    before = nrow(df)
    dropmissing!(df, cols)
    after = nrow(df)
    if before > after
        @warn "dropping rows: $before -> $after" colHasMiss
    end
end

dfs = readdf.(args.infiles);
@info "infile(s):" size.(dfs)
df = joinleft(dfs, args.join);
@info "after leftjoins:" size(df)
df = dropmissing_cols(df)

# the first class column found in training data is the class used for training.
class = intersect(args.class, names(df))[1]
@info "class for training: \"$class\""

if args.test !== nothing
    if isfile(args.test)
        df_test = readdf(args.test);
        @info "test:" size(df_test)
        df_test = joinleft([df_test, dfs[2:end]...], args.join);
        @info "after leftjoins:" size(df_test)
        df_test = dropmissing_cols(df_test)
        # Any columns in df_test not found in df will not be features, and will 
        # be written out unaffected. Extra columns in train should only be the 
        # -c/--class or non-number since we currently ignore them.
        extraInTrain = setdiff(names(df), names(df_test))
        extraInTrainTypes = df[!, extraInTrain] |> eachcol .|> eltype
        # we are only concerned about potential features
        extraInTrain = extraInTrain[extraInTrainTypes .<: Real]
        # class can be missing from test data, we might not know the truth 
        # about something we are predicting on.
        extraInTrain = extraInTrain[extraInTrain .!= class]
        @assert length(extraInTrain) == 0 "Columns missing in -t/--test: $extraInTrain"
    elseif args.test ∈ names(df)
        df_test = df[args.test .== 1]
        df = df[args.test .!= 1]
        select!(df_test, Not(args.test))
        select!(df, Not(args.test))
    else
        try nTest = parse(Int, args.test)
        catch
            try nTest = parse(Float64, args.test) * nrow(df)
            catch
                error("Not understood: -t/--test " + args.test)
            end
        end
        testIdx = shuffle([ones(Int, nTest); zeros(Int, nrow(df) - nTest)])
        df_test = df[testIdx]
        df = df[.!testIdx]
    end
end

"Flatten per sequence."
function substitute(sequence::Vector{Char}, alphabet::Vector{Char}, sub::Matrix)
    dropdims(vcat((sub[:, alphabet .== l] for l in sequence)...); dims=2)
end

# perform substitution
if args.substitute !== nothing
    subCol, subPath = args.substitute
    if !endswith(subPath, ".tsv") subPath *= ".tsv" end
    if !isfile(subPath) subPath = "$ROOT/data/NCBI/" * subPath end
    @assert isfile(subPath) "File not found: $(args.substitute[2])"
    subMat = readdf(subPath)
    # gap is written as both * and -
    if "-" ∉ names(subMat) && "*" ∈ names(subMat) subMat[!, "-"] .= subMat[!, "*"] end
    if "*" ∉ names(subMat) && "-" ∈ names(subMat) subMat[!, "*"] .= subMat[!, "-"] end
    # find single string column
    stringColInd = findall(eltype.(eachcol(subMat)) .<: AbstractString) |> only
    outAlphabet = subMat[!, stringColInd]
    subMat = subMat[!, Not(stringColInd)]
    inAlphabet = names(subMat) .|> only
    subMat = Matrix(subMat)
    
    function substitute(df::DataFrame)
        # make dataframe mapping from unique sequences to their substitutions 
        sequences = unique(df[!, subCol])
        seqlen = length.(sequences) |> unique |> only
        decimals = seqlen |> string |> length
        newSubCols = subCol*"_" .* vcat([outAlphabet .* lpad(i, decimals, '0') for i in 1:seqlen]...)
        subs = hcat(substitute.(collect.(sequences), Ref(inAlphabet), Ref(subMat))...)
        df_sub = DataFrame(subs', newSubCols)
        df_sub[!, subCol] = sequences
        
        # annotate the data with the substitutions. Inner- instead of leftjoin 
        # since it is guranteed that all sequences will be found, since they 
        # are taken from df and this will mean we get Float64 instead of 
        # Float64? (allows missing).
        df = innerjoin(df, df_sub; on=subCol)
        select!(df, Not(subCol))
        @info "after substitutions:" size(df)
        df
    end
    df = substitute(df);
    if args.test !== nothing df_test = substitute(df_test); end
end

# drop trivial features
trivial = @chain df Matrix all(_ .== _[[1], :]; dims=1) vec
if any(trivial)
    @info "Dropping trivial" names(df)[trivial]
    df = df[!, .!trivial];
end

X = df[!, Not(args.join)]
numCols = eltype.(eachcol(X)) .<: Union{Missing,Number}
if !all(numCols)
    @info "Dropping non-number columns: $(names(X)[.!numCols])"
    X = X[!, numCols];
end;
if args.test !== nothing
    Xnew = df_test[!, unique([names(X); intersect(args.class, names(df_test))])];
else
    Xnew = DataFrame(zeros(0, size(X,2)), names(X))
end

dropmissing_verbose!(X)
dropmissing_verbose!(Xnew)

y = X[!, class]
X = X[!, Not(class)]
ynew = Xnew[!, intersect(args.class, names(Xnew))]
Xnew = Xnew[!, names(X)];

function Base.parse(x::AbstractString)
    try parse(Int, x)
    catch
        try parse(Float64, x)
        catch; x end
    end
end

getpred(clf, Xnew) = predict_proba(clf, Matrix(Xnew))[:, 2]


function evalpred(clf, Xnew, ynew::Vector)
    pred = getpred(clf, Xnew)
    
    tru = ynew
    if length(unique(tru)) == 2
        println("AUC = ", AUC(pred, tru))
    else
        println("cor = ", cor(pred, tru))
        println("AUC = ", AUC(pred, tru .> mean(tru)))
    end
end
evalpred(clf, Xnew, ynew::DataFrame) = evalpred(getpred(clf, Xnew), ynew)
function evalpred(pred, ynew::DataFrame)
    for class_test in names(ynew)
        println(class_test)
        tru = ynew[!, class_test]
        if length(unique(tru)) == 2
            println("AUC = ", AUC(pred, tru))
        else
            println("cor = ", cor(pred, tru))
            println("AUC = ", AUC(pred, tru .> mean(tru)))
        end
    end
end

model = Chain(Dense(size(X,2) => 16, tanh), Dense(16 => 1, σ), vec)

_X = Matrix{Float32}(X)' |> collect
optim = Flux.setup(Flux.Adam(), model);  # will store optimiser momentum, etc.

for epoch in 1:2_000
    loss, grads = Flux.withgradient(model) do m
        preds = m(_X)
        sum(preds .^ 2 .* (1 .- 2 .* y))
    end
    Flux.update!(optim, model, grads[1])
    epoch % 100 == 0 && println(loss)
end

preds = model(Matrix{Float32}(Xnew)')
AUC(preds, ynew.reaction)
CSV.write("hejsa.tsv", DataFrame(reaction=ynew.reaction, pred=preds); delim='\t')



global clf
if args.train !== nothing
    model_args = [Symbol(k)=>parse(v) for (k,v) in split.(args.model, '=')]
    clf = RandomForestClassifier(; model_args...)
    
    @time fit!(clf, Matrix(X), y)
    JLD.save(args.train, "model", clf)
else
    # clf = JLD.load(args.eval, "model")
    # in case a model is saved with another key
    clf = JLD.load(args.eval) |> values |> only
    evalpred(clf, X, y)
end

nrow(Xnew) == 0 || evalpred(clf, Xnew, ynew)




