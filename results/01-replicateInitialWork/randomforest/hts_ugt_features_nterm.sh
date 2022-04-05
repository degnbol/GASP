#!/usr/bin/env zsh
mlr --tsv --from hts_ugt_features.tsv cut -x -f seq then join -j enzyme -f muscle_nterm.tsv > hts_ugt_features_nterm.tsv
