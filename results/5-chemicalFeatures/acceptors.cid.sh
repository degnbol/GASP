#!/usr/bin/env zsh
ROOT=`git root`
mlr --tsv uniq -f cid $ROOT/results/*-validateAcceptors/reactions.tsv \
    $ROOT/data/acceptorsOfInterest/acceptors-of-interest.tsv | sed 1d |
    cat - $ROOT/data/inhouse/david/negatives.cid |
    grep -v '^$' | sort -u > acceptors.cid
