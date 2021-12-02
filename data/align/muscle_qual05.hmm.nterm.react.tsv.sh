#!/usr/bin/env zsh
# filter consensus alignment by enzymes that are found among the reactivity data.
mlr --tsv --from ../reactions.tsv uniq -f enzyme then join -j enzyme -f muscle_qual05.hmm.nterm.tsv > muscle_qual05.hmm.nterm.react.tsv
