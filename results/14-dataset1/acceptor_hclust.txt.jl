#!/usr/bin/env julia
# hierarchical clustering of chemicals for visuals.
using Clustering
using DataFrames, CSV
using Statistics
using Distances

df = CSV.read("dataset1.tsv", DataFrame)

df_mat = unstack(df, :acceptor, :enzyme, :rate)
mat = Matrix(df_mat[!, 2:end])

# replace missing values with average rate per enzyme
avg_enz_rates = mat |> eachcol .|> skipmissing .|> mean
for (col, avg_rate) in enumerate(avg_enz_rates)
    mat[ismissing.(mat[:, col]), col] .= avg_rate
end

dists = pairwise(Euclidean(), mat; dims=1)

ordering = hclust(dists; linkage=:ward, branchorder=:optimal).order
ordered = df_mat.acceptor[ordering]

open("acceptor_hclust.txt", "w") do io
    println(io, join(ordered, '\n'))
end

