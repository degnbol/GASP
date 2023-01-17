#!/usr/bin/env julia
using CSV, DataFrames
using Chain: @chain
using Random
using Statistics
using LinearAlgebra
using DecisionTree
using ThreadPools
using Printf
ROOT = readchomp(`git root`)
include("$ROOT/src/utilities/glob.jl")
include("$ROOT/src/utilities/AUC.jl")

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
function topP_metric(pred::Vector{Float64}, class::Union{BitVector,Vector{Bool}})::Float64
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


missmax(x) = x[.!ismissing.(x)] |> maximum


### READ

df = @chain "$ROOT/results/*-features/train.tsv" glob CSV.read(DataFrame) unique
df_seq_feat = @chain "$ROOT/results/*-features/blosum62Amb.tsv.gz" glob CSV.read(DataFrame)
df_chem_feat = @chain "$ROOT/results/*-chemicalFeatures/acceptors2_features.tsv" glob CSV.read(DataFrame)

df.reaction = df.reaction .|> Bool 

leftjoin!(df, df_seq_feat; on=:enzyme)
leftjoin!(df, df_chem_feat; on=:cid)

before = nrow(df)
dropmissing!(df)
println("dropmissing: $before -> ", nrow(df))

nSample = 2000
println("sample $nSample datapoints from generated negatives.")
isGenNeg = startswith.(df.source, "negatives")
df_negSample = df[isGenNeg, :]
df_negSample = df_negSample[shuffle([trues(nSample); falses(nrow(df_negSample)-nSample)]), :]
df = vcat(df[.!isGenNeg, :], df_negSample)

feat_names = setdiff(names(df), ["reaction", "enzyme", "cid", "source"])
realcols = eltype.(eachcol(df[!, feat_names])) .<: Real
if !all(realcols)
    @info "Dropping non-number columns: $(feat_names[.!realcols])"
    feat_names = feat_names[realcols]
end;

df_feats = DataFrame(name=feat_names)
df_feats.unique = df[!, feat_names] |> eachcol .|> unique .|> length
df_feats.var = df[!, feat_names] |> eachcol .|> var
trivial = df_feats.name[df_feats.unique .== 1]
if length(trivial) > 0
    @info "Dropping trivial" trivial
    df = df[!, Not(trivial)];
    setdiff!(feat_names, trivial)
end;

# correlation between features
feat_cors = @chain df[!, feat_names] Matrix cor
feat_cors[diagind(feat_cors)] .= 0.; # ignore autocorr
feat_maxcors = maximum(feat_cors; dims=2) |> vec
println("Perfectly correlated features: ", sum(feat_maxcors .== 1))
feat_maxcor_names = [feat_names[idx] for idx in eachrow(feat_cors .== feat_maxcors)]
df_feats_cor = DataFrame(name=feat_names, maxcor=feat_maxcors, maxcor_names=join.(feat_maxcor_names, ' '))
leftjoin!(df_feats, df_feats_cor; on=:name)

df_feats.iteration .= 0

# only consider one of a set of perfectly correlated features.
ignore = String[]
for row in eachrow(df_feats_cor[df_feats_cor.maxcor .== 1, [:name, :maxcor_names]])
    row.name in ignore || append!(ignore, split(row.maxcor_names, ' '))
end
unique!(ignore)
println("Ignoring $(length(ignore)) perfectly correlated features.")
setdiff!(feat_names, ignore)

isSeqFeat = startswith.(feat_names, "seq_")

for (colname, consider_idx, n_trees) in [(:cid, .!isSeqFeat, 100), (:enzyme, isSeqFeat, 20)]
    consider = feat_names[consider_idx]

    uCol = unique(df[!, colname])
    nConsider = length(consider)

    X, y = Matrix{Float64}(df[!, consider]), df[!, :reaction]

    for it in 1:nConsider-1
        sample = shuffle(uCol)[1:floor(Int, length(uCol) / 5)]
        testidx = df[!, colname] .∈ Ref(sample)

        Xtrain, ytrain = X[.!testidx, :], y[.!testidx]
        Xtest, ytest = X[testidx, :], y[testidx]

        @time metrics = tmap(1:size(X,2)) do i
            clf = RandomForestClassifier(n_trees=n_trees)
            fit!(clf, view(Xtrain, :, Not(i)), ytrain)
            pred = predict_proba(clf, view(Xtest, :, Not(i)))[:, 2]
            topP_metric(pred, ytest)
        end

        leftjoin!(df_feats, DataFrame("name"=>consider, @sprintf("metric_%s_%03d", colname, it)=>metrics); on=:name)

        best = argmax(metrics)
        df_feats.iteration[df_feats.name .== consider[best]] .= it
        println("Selected: ", consider[best])
        flush(stdout)
        deleteat!(consider, best)
        X = X[:, Not(best)]
    end

    # the winner:
    df_feats.iteration[df_feats.name .== only(consider)] .= nConsider
    
    CSV.write("feature_selection.tsv.gz", df_feats; delim='\t', compress=true)
    
    max_metrics = df_feats[!, startswith.(names(df_feats), "metric_$colname")] |> eachcol .|> missmax
    max_metric, max_it = findmax(max_metrics)
    max_selected = [only(df_feats.name[df_feats.iteration .== i]) for i in 1:max_it]
    println(@sprintf("Max metric %.4f reached with removal of: ", max_metric), join(max_selected, ", "))
end

