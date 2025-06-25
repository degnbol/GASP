#!/usr/bin/env julia
using Pkg
Pkg.instantiate()

# E3FP
# make sure to run this with the correct env active.
# This is to build the PyCall with an env that works for E3FP python package.
ENV["PYTHON"] = readchomp(`which python`)
Pkg.build("PyCall")

