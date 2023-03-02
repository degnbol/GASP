#!/usr/bin/env zsh
ROOT=`git root`
FEATURES=`\ls -d $ROOT/results/*features`
ids="enzyme acceptor source cid"
chemFeatures=$(gzcat $FEATURES/traintest_match.tsv | head -n1 | tr '\t' '\n' | grep -v 'seq_' | tr '\n' ' ')
features=$(cat $FEATURES/selection/species_select.txt | tr '\n' ' ')
gzcat $FEATURES/traintest_match.tsv.gz | randomforest.py -t 1 -i ${=ids} -c reaction -e rate -n 10000 -F ${=chemFeatures} ${=features} > pred_onehot-speciesSeqSelect.tsv
gzcat $FEATURES/traintest_blosum62.tsv.gz | randomforest.py -t 1 -i ${=ids} -c reaction -e rate -n 10000 -F ${=chemFeatures} ${=features} > pred_blosum-speciesSeqSelect.tsv
