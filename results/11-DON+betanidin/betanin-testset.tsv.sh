#!/usr/bin/env zsh
mlr -t template -f enzyme,betaninProd + put '$cid = 135449343' \
    `git root`/data/inhouse/david/betaninProd.tsv ../07-align/muscle_qual.hmm.nterm.tsv.gz > betanin-testset.tsv
