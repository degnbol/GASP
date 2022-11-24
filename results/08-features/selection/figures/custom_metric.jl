#!/usr/bin/env julia
using Gadfly, Cairo, Fontconfig
using DataFrames
using Random

n = 10
npos = 4
nneg = n-npos
nTP = 2
nFP = npos-nTP

df = DataFrame(pred=rand(n))
sort!(df, :pred; rev=true)

df[!, :weight] .= [npos:-1:1; zeros(n-npos)]
df[!, :weight] ./= sum(df.weight)

pos = shuffle([fill(true, nTP); fill(false, npos-nTP)])
neg = shuffle([fill(true, nFP); fill(false, nneg-nFP)])
df[!, :class] .= [pos; neg]

p = plot(df, x=:pred, y=:weight, color=:class, Geom.point, Geom.hair,
         Scale.color_discrete_manual("lightgray", "lime"));

p |> PDF("custom_metric.pdf", 10cm, 8cm)

