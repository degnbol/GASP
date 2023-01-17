#!/usr/bin/env zsh
ROOT=`git root`
mlr -t --from experimentYield.tsv join -j enzyme -f $ROOT/results/*-align/muscle_qual.hmm.nterm.tsv.gz + join -j cid --lk smiles -f /Users/cdmadsen/GT/results/05-chemicalFeatures/acceptors2_features.tsv + put '$Enzymes = gsub($seq, "-", "")' > forSubFinder.tsv
mlr --t2c --from forSubFinder.tsv cut -x -f Yield + rename smiles,Metabolites + split -n 500 --prefix=forSubFinder
