#!/usr/bin/env zsh
ids="enzyme acceptor source cid"
zcat traintest_chem_blosum62Amb.tsv.gz |
    randomforest.py -t 1 -i ${=ids} -c reaction -e rate -n 10000 |
    gzip > pred_chem_blosum62Amb.tsv.gz

