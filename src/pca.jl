#!/usr/bin/env julia
using PyCall
@pyimport numpy as np
@pyimport scipy.sparse as sp
using LinearAlgebra
using SparseArrays
using Arpack

fname = expanduser("~/biosustain/gt/acceptors/acceptors.fpz")
# taken from python  package e3fp.fingerprint.db.FingerprintDatabase
fpdb = Dict(np.load(fname, allow_pickle=true).items())

X = SparseMatrixCSC{Bool,Int}(fpdb["shape"][2], fpdb["shape"][1], fpdb["indptr"], fpdb["indices"], fpdb["data"])

svds(X'; nsv=3, tol=1e-1, maxiter=2)

