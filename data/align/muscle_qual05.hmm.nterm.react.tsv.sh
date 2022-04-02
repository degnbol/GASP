#!/usr/bin/env zsh
# filter consensus alignment by enzymes that are found among the reactivity data.
INFILE=`git root`/data/reactions/reactions.tsv
mlr --tsv --from $INFILE uniq -f enzyme then join -j enzyme -f muscle_qual05.hmm.nterm.tsv > muscle_qual05.hmm.nterm.react.tsv
