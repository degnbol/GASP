#!/usr/bin/env zsh
ROOT=`git root`
mlr -t --from $ROOT/data/reactions/reactions.tsv filter '$source == "pTMH" || $source == "HTS_UGT"' +\
    join -j enzyme -f $ROOT/results/*-align/muscle_qual.hmm.nterm.tsv.gz + uniq -a > dataset1.tsv
