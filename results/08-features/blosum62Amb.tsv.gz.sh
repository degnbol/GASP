#!/usr/bin/env zsh
# simply encode all sequences with blosum62 (with ambiguous codes).
# All features were kept, but one could use -k flag to ensure.
ROOT=`git root`
gunzip -c $ROOT/results/*-align/muscle_qual.hmm.nterm.tsv.gz |
    encode_features.py -i enzyme --aa seq --aa-encoding $ROOT/data/NCBI/blosum62Amb.tsv |
    gzip -c > blosum62Amb.tsv.gz
