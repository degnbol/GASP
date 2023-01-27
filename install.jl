#!/usr/bin/env julia
# use conda env 
using Pkg

Pkg.add(["ArgParse", "DelimitedFiles"])
# for MDS
Pkg.add("MultivariateStats")
# for traintest.jl
Pkg.add(["DataFrames", "CSV", "Random", "MLJ", "Glob", "DecisionTree", "Chain", "ThreadPools"])
Pkg.add(url="https://github.com/diegozea/ROC.jl")

Pkg.add(["Conda", "PyCall"])
Pkg.build("PyCall")
using Conda
Conda.add("e3fp")


