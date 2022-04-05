#!/usr/bin/env zsh
# filter consensus alignment by enzymes that are found among the reactivity data.
ROOT=`git root`
INFILE=`ls $ROOT/results/*validateAcceptors/reactions.tsv`
mlr --tsv --from $INFILE uniq -f enzyme then join -j enzyme -f muscle_qual05.hmm.nterm.tsv > muscle_qual05.hmm.nterm.react.tsv
