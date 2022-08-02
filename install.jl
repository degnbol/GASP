#!/usr/bin/env julia
# use conda env 
using Pkg

Pkg.add(["MultivariateStats", "ArgParse", "DelimitedFiles"])

ENV["PYTHON"] = expanduser("~/miniconda3/envs/gt/bin/python")
Pkg.add("PyCall")
# if PyCall was already installed for some reason then explicitly calling build 
# will be necessary after setting the PYTHON env variable.
Pkg.build("PyCall")
# Now PyCall should be using the python env where we installed things like 
# E3FP. Julia may have to be restarted for E3FP to be found with 
# pyimport("E3FP")

