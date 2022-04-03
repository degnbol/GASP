#!/usr/bin/env zsh
ids="enzyme acceptor source cid"
zcat traintest_chem_matchAmb.tsv.gz |
    randomforest.py -t 1 -i ${=ids} -c reaction -e rate -n 10000 > pred_chem_matchAmb.tsv

