mlr --tsv --from gtpred_features.tsv cut -x -f seq then join -f gtpred.muscle.hmm.tsv -j enzyme > gtpred_features_hmm.tsv
