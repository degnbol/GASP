#!/usr/bin/env zsh
mlr --tsv --from gtpred_features.tsv cut -x -f seq then join -j enzyme -f muscle_nterm.tsv > gtpred_features_nterm.tsv
