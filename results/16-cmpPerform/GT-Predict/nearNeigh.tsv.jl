#!/usr/bin/env julia
using CSV, DataFrames
ROOT = readchomp(`git root`)
include("$ROOT/src/utilities/glob.jl")

## READ

df_water = CSV.read("water.tsv.gz", DataFrame)
df_gtpred = CSV.read("$ROOT/data/reactions/gtpred_reactions_enz.tsv", DataFrame)
df_react = CSV.read(glob("$ROOT/results/*generateNegatives/reactions.tsv"), DataFrame)

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

# For each query enzyme, we use values from either nearest neighbor or second 
# nearest neighbor when there is no entry from nearest neighbor and a given
# acceptor.












