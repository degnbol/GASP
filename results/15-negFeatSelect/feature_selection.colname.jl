#!/usr/bin/env julia
using CSV, DataFrames

df = CSV.read("feature_selection.tsv.gz", DataFrame)
df = df[.!startswith.(df.name, "seq_"), :]
@info "$(nrow(df)) features in total"

missmax(x) = x[.!ismissing.(x)] |> maximum

selected = Dict{String,Vector{String}}()

for metric_name in ["topP", "AUC"]
    max_metrics = df[!, startswith.(names(df), metric_name)] |> eachcol .|> missmax
    max_metric, max_it = findmax(max_metrics)
    max_deselected = [only(df.name[df.iteration .== i]) for i in 1:max_it]
    @info "$(length(max_deselected)) features deselected for max $metric_name"
    # any perfectly correlated features would have been ignored in the analysis 
    # and would yield exactly the same performance if the deselected twin had 
    # been ignored instead so we also deselect those.
    twins = vcat(split.(df[df.name .∈ Ref(max_deselected) .&& df.maxcor .== 1, :maxcor_names], ' ')...)
    @info "deselecting a further $(length(twins)) features: $(join(twins, ','))"
    append!(max_deselected, twins)
    selected[metric_name] = setdiff(df.name, max_deselected)
    open("feature_selection-$metric_name.colname", "w") do io
        println(io, join(selected[metric_name], '\n'))
    end
end

