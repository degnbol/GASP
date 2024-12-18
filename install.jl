#!/usr/bin/env julia
# use conda env 
using Pkg

Pkg.add(["ArgParse", "DelimitedFiles"])
# for MDS
Pkg.add("MultivariateStats")
# for traintest.jl
Pkg.add(["DataFrames", "CSV", "Random", "MLJ", "Glob", "DecisionTree", "Chain", "ThreadPools"])
Pkg.add(url="https://github.com/diegozea/ROC.jl")

# E3FP
Pkg.add("PyCall")
# make sure to run `conda activate GT` before this.
# This is to build the PyCall with an env that works for E3FP python package.
ENV["PYTHON"] = readchomp(`which python`)
Pkg.build("PyCall")


