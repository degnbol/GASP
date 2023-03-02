#!/usr/bin/env zsh
ROOT=`git root`
FEATURES=`\ls -d $ROOT/results/*features`
ids="enzyme acceptor source cid"
seqFeatures=$(gzcat $FEATURES/traintest_match.tsv | head -n1 | tr '\t' '\n' | grep 'seq_' | tr '\n' ' ')
features=$(cat $FEATURES/selection/selected_ratethres.txt | tr '\n' ' ')
gzcat $FEATURES/traintest_match.tsv | randomforest.py -t 1 -i ${=ids} -c reaction -e rate -n 10000 -F ${=seqFeatures} ${=features} > pred_onehot-chemChemSelect.tsv
gzcat $FEATURES/traintest_blosum62.tsv | randomforest.py -t 1 -i ${=ids} -c reaction -e rate -n 10000 -F ${=seqFeatures} ${=features} > pred_blosum-chemChemSelect.tsv
