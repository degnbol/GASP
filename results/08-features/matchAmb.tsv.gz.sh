#!/usr/bin/env zsh
# simply encode all sequences with match (with ambiguous codes).
# Less than 2% of features never vary and were discarded, this can be avoided with flag -k.
ROOT=`git root`
cat $ROOT/results/07-align/muscle_qual.hmm.nterm.tsv |
    encode_features.py -i enzyme --aa seq --aa-encoding $ROOT/data/NCBI/matchAmb.tsv |
    gzip -c > matchAmb.tsv.gz
