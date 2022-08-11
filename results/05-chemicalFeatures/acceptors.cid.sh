#!/usr/bin/env zsh
ROOT=`git root`
mlr -t uniq -f cid $ROOT/results/*-generateNegatives/reactions.tsv |
    sed 1d | sort -u > acceptors.cid
mlr -t uniq -f cid $ROOT/data/acceptorsOfInterest/acceptors-of-interest.tsv |
    sed 1d >> acceptors.cid

