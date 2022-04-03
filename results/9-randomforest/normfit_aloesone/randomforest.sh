#!/usr/bin/env zsh
ids="enzyme acceptor source cid"
gzcat `git root`/results/*features/traintest_matchAmb.tsv.gz | randomforest.py -t 1 -i ${=ids} -c reaction -e rate -n 10000 > pred_onehot.tsv

