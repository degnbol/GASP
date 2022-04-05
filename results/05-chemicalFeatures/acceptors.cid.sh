#!/usr/bin/env zsh
ROOT=`git root`
mlr --tsv uniq -f cid $ROOT/results/*-validateAcceptors/acceptors.tsv |
    sed 1d | sort -u > acceptors.cid
