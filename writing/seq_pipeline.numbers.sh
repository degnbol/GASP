#!/usr/bin/env zsh
# print some numbers for the seq pipeline figure.
sed 1d `git root`/results/*-align/muscle_qual.hmm.nterm.tsv | wc -l | xargs | tr -d '\n'
echo " consensus alignments"
