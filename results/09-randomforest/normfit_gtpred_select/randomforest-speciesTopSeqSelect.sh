#!/usr/bin/env zsh
ROOT=`git root`
FEATURES=`\ls -d $ROOT/results/*features`
ids="enzyme acceptor source cid"
chemFeatures=$(gzcat $FEATURES/traintest_match.tsv | head -n1 | tr '\t' '\n' | grep -v 'seq_' | tr '\n' ' ')
features=$(cat $FEATURES/selection/speciesTop_select.txt | tr '\n' ' ')
gzcat $FEATURES/traintest_match.tsv | randomforest.py -t 1 -i ${=ids} -c reaction -e rate -n 10000 -F ${=chemFeatures} ${=features} > pred_onehot-speciesTopSeqSelect.tsv
gzcat $FEATURES/traintest_blosum62.tsv | randomforest.py -t 1 -i ${=ids} -c reaction -e rate -n 10000 -F ${=chemFeatures} ${=features} > pred_blosum-speciesTopSeqSelect.tsv
