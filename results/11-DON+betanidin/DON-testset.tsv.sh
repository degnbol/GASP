#!/usr/bin/env zsh
mlr -t --from ../07-align/muscle_qual.hmm.nterm.tsv.gz uniq -f enzyme + \
    put '$cid = 40024' > DON-testset.tsv
