#!/usr/bin/env zsh
ROOT=`git root`
FEATURES=`\ls -d $ROOT/results/*features`
ids="enzyme acceptor source cid"
# use features selected for distinguishing chemical groups AND for distinguising between species
features=$(cat $FEATURES/selection/species_select.txt | tr '\n' ' ')
gzcat $FEATURES/traintest_match.tsv | randomforest.py -t 1 -i ${=ids} -c reaction -e rate -n 10000 -F ${=features} > pred_onehot.tsv
gzcat $FEATURES/traintest_blosum62.tsv | randomforest.py -t 1 -i ${=ids} -c reaction -e rate -n 10000 -F ${=features} > pred_blosum.tsv
