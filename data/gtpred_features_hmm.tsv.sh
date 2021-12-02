#!/usr/bin/env zsh
mlr --tsv --from gtpred_features.tsv cut -x -f seq then join -f hmm_alignment/gtpred.muscle.hmm.tsv -j enzyme > gtpred_features_hmm.tsv
