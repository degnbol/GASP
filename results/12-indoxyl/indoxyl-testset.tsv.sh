#!/usr/bin/env zsh
ROOT=`git root`
mlr -t --from $ROOT/results/*-align/muscle_qual.hmm.nterm.tsv.gz uniq -f enzyme + \
    put '$cid = 50591' > indoxyl-testset.tsv
