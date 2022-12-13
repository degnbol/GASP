#!/usr/bin/env zsh
ROOT=`git root`
train=$ROOT/results/08-features/train.tsv 
mlr -t --from $train uniq -f cid + join -j cid -f experimentYield.tsv > testsetNovelSeqs.tsv

