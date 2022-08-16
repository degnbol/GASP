#!/usr/bin/env julia
ROOT = readchomp(`git root`)
include("$ROOT/src/utilities/glob.jl")
using DataFrames, CSV
using Statistics
using Chain: @chain
readtsv(path::AbstractString) = CSV.read(path, DataFrame)

seqs = readtsv("$ROOT/results/07-align/muscle_qual.hmm.nterm.tsv")

testset = @chain "../pred_all-aucEnz.tsv" readtsv leftjoin(seqs; on=:enzyme)

trainset = @chain "$ROOT/results/*-features/train.tsv" begin
    glob
    only
    readtsv
    groupby(:enzyme)
    combine(nrow => :N)
    leftjoin(seqs; on=:enzyme)
    dropmissing
end

testseqs = hcat(collect.(testset[!, :seq])...)
trainseqs = hcat(collect.(trainset[!, :seq])...)
# reshape to (seqLen, 1, nTrain) so it can compare to (seqLen, nTest, [implicit 1])
trainseqs = @chain trainseqs cat(; dims=3) permutedims((1, 3, 2))
nMatches = @chain testseqs .== trainseqs sum(; dims=1) dropdims(;dims=1)

seqLen = size(testseqs, 1)
distances = (seqLen .- nMatches) ./ seqLen

testset[!, :minDist] .= minimum(distances; dims=2)
testset[!, :meanDist] .= mean(distances; dims=2)
testset[!, :medianDist] .= median(distances; dims=2)

rename!(testset,
        :Yield_thres_50_pred_AUC => :auc,
        :Yield_thres_50_pred_J_thres_acc => :acc,
        :Yield_pred_corr => :cor)

testset2 = dropmissing(testset)

Statistics.cor(X::DataFrame, Y::DataFrame) = cor(Matrix(X), Matrix(Y))
Statistics.cor(X::Vector, Y::DataFrame) = cor(X, Matrix(Y))

distNames = [:minDist, :meanDist, :medianDist]
cor(testset[!, [:acc, :cor]], testset[!, distNames])
cor(testset2[!, :auc], testset2[!, distNames])

using Plots

scatter(testset[!, :minDist], testset[!, :acc]; ylab="acc",
        xlim=[0,1], ylim=[0,1], xlab="dist", label="min", legend=:outertopright)
scatter!(testset[!, :meanDist], testset[!, :acc]; label="mean")
scatter!(testset[!, :medianDist], testset[!, :acc]; label="median")
savefig("acc.pdf")

scatter(testset[!, :minDist], testset[!, :cor]; ylab="cor",
        xlim=[0,1], ylim=[0,1], xlab="dist", label="min", legend=:outertopright)
scatter!(testset[!, :meanDist], testset[!, :cor]; label="mean")
scatter!(testset[!, :medianDist], testset[!, :cor]; label="median")
savefig("cor.pdf")

scatter(testset2[!, :minDist], testset2[!, :auc]; ylab="auc",
        xlim=[0,1], ylim=[0,1], xlab="dist", label="min", legend=:outertopright)
scatter!(testset2[!, :meanDist], testset2[!, :auc]; label="mean")
scatter!(testset2[!, :medianDist], testset2[!, :auc]; label="median")
savefig("auc.pdf")

