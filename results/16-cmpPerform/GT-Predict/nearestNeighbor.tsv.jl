#!/usr/bin/env julia
using CSV, DataFrames
using StatsBase: ordinalrank, mean
using Chain: @chain
using Printf
ROOT = readchomp(`git root`)
include("$ROOT/src/utilities/glob.jl")
include("$ROOT/src/utilities/AUC.jl")

## READ

df_water = CSV.read("water.tsv.gz", DataFrame)
# the following file is the same as the raw GT-Predict 
# enzyme_interaction_data.txt except for some formatting and rename UGT -> At_
df_gtpred = CSV.read("$ROOT/data/reactions/gtpred_reactions_enz.tsv", DataFrame)
df_react = CSV.read(glob("$ROOT/results/*generateNegatives/reactions.tsv"), DataFrame)
df_raw2cid = CSV.read("$ROOT/data/reactions/rawAcceptor_cid_title.tsv", DataFrame)

## PREPROCESS

df_react = df_react[df_react.reaction .!= 0.5, :]
df_react.reaction = df_react.reaction .|> Bool

df_gtpred = stack(df_gtpred, Not(:acceptor), variable_name=:enzyme)
# GT-Predict values as per
# GTPredict+code/PredictEnzymeInteraction/Readme.txt
# and as per the source code in
# GTPredict+code/PredictEnzymeInteraction/matlab\ code/PredictEnzymeInteraction.m
# - missing = negative
# - 0       = negative
# - 1       = positive
# - 2       = unknown
# - 3       = missing
# missing entries are set to zero in the code since the table is read with 
# dlmread https://au.mathworks.com/help/matlab/ref/dlmread.html in
# GTPredict+code/PredictEnzymeInteraction/matlab\ 
# code/read_enzyme_interaction_file.m
df_gtpred.value[ismissing.(df_gtpred.value)] .= 0
df_gtpred = df_gtpred[df_gtpred.value .< 2, :]
df_gtpred.value = df_gtpred.value .|> Bool

## NEAREST NEIGHBOR

# For each query enzyme, we use values from either nearest neighbor or second 
# nearest neighbor when there is no entry from nearest neighbor and a given
# acceptor.

df_water = transform(groupby(df_water, :Query), :Score => (x->ordinalrank(x,rev=true)) => :rank)
df_water = df_water[df_water.rank .<= 2, :]
# Rename to match column for joining
df_water.Target = replace.(df_water.Target, "UGT" => "At_")

## JOINS

df = innerjoin(df_water, df_gtpred, on= :Target => :enzyme)

# nearest neighbor, or second nearest if no first nearest but never third or 
# higher.
df = combine(groupby(df, [:Query, :acceptor])) do sdf
    sort(sdf, :Score, rev=true) |> first
end

rename!(df, :Query => :enzyme, :Target => :Neighbor, :acceptor => :raw, :value => :GTpred)

df_chem = unique(df[!, [:raw]])
leftjoin!(df_chem, df_raw2cid, on=:raw)
@info "Unmapped chems:" df_chem[ismissing.(df_chem.cid), :raw]

leftjoin!(df, df_chem; on=:raw)
df = innerjoin(df, df_react; on=[:enzyme, :cid], matchmissing=:notequal)
CSV.write("nearestNeighbor.tsv", df; delim='\t')

acc(ŷ, y) = mean(ŷ .== y)

@printf "Incl identical: acc = %.3f, AUC = %.3f\n" acc(df.GTpred, df.reaction) AUC(df.GTpred, df.reaction)
@printf "Excl identical: acc = %.3f, AUC = %.3f\n" [f(eachcol(df[df.enzyme .!= df.Neighbor, [:GTpred, :reaction]])...) for f in [acc, AUC]]...
@printf "Only identical: acc = %.3f, AUC = %.3f\n" [f(eachcol(df[df.enzyme .== df.Neighbor, [:GTpred, :reaction]])...) for f in [acc, AUC]]...

