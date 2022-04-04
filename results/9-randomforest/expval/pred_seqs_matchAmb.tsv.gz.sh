#!/usr/bin/env zsh
ids="enzyme acceptor source cid"
gzcat traintest_seqs_matchAmb.tsv.gz |
    randomforest.py -t 1 -i ${=ids} -c reaction -e rate -n 10000 |
    gzip > pred_seqs_matchAmb.tsv.gz

