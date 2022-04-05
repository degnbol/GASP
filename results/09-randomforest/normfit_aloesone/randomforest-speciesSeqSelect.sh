#!/usr/bin/env zsh
features=`git root`/results/*features
ids="enzyme acceptor source cid"
chemFeatures=$(gzcat $features/traintest_match.tsv | head -n1 | tr '\t' '\n' | grep -v 'seq_' | tr '\n' ' ')
features=$(cat $features/selection/species_select.txt | tr '\n' ' ')
gzcat $features/traintest_matchAmb.tsv.gz | randomforest.py -t 1 -i ${=ids} -c reaction -e rate -n 10000 -F ${=chemFeatures} ${=features} > pred_onehot-speciesSeqSelect.tsv

