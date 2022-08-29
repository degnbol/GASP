#!/usr/bin/env julia
using ROC # https://github.com/diegozea/ROC.jl
AUC(ŷ, y) = ROC.roc(ŷ, y) |> ROC.AUC
