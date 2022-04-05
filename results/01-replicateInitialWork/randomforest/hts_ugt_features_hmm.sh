mlr --tsv --from hts_ugt_features.tsv cut -x -f seq then join -f hts_ugt-gtpred.muscle.hmm.tsv -j enzyme > hts_ugt_features_hmm.tsv
