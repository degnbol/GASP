#!/usr/bin/env zsh
ROOT=`git root`
seqs=$ROOT/results/07-align/muscle_qual.hmm.nterm.tsv
train=$ROOT/results/08-features/train.tsv 

mlr -t --from experimentYield.tsv join -j enzyme -f $seqs > testset.tmp
mlr -t --from $train uniq -f enzyme + join -j enzyme -f $seqs + join -j seq -f testset.tmp --np --ul > testsetDist.tsv
rm testset.tmp

