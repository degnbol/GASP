#!/usr/bin/env zsh
FEATURES=`git root`/results/*features
ids="enzyme acceptor source cid"
gzcat $FEATURES/traintest_match.tsv.gz | randomforest.py -t 1 -i ${=ids} -c reaction -e rate -n 10000 > pred_onehot.tsv
gzcat $FEATURES/traintest_blosum62.tsv.gz | randomforest.py -t 1 -i ${=ids} -c reaction -e rate -n 10000 > pred_blosum.tsv
