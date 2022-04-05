#!/usr/bin/env zsh
FEATURES=`git root`/results/*features
ids="enzyme acceptor source cid"
features=$(tr '\n' ' ' < $FEATURES/selection/selected_ratethres.txt)
randomforest.py -t 1 -i ${=ids} -c reaction -e rate -n 10000 -F ${=features} < $FEATURES/test_negatives_matchAmb.tsv > pred_onehot.tsv
randomforest.py -t 1 -i ${=ids} -c reaction -e rate -n 10000 -F ${=features} < $FEATURES/test_negatives_blosum62Amb.tsv > pred_blosum.tsv
mlr --tsv -I --from pred_onehot.tsv cut -f cid,enzyme,pred
mlr --tsv -I --from pred_blosum.tsv cut -f cid,enzyme,pred
table.py pred_onehot.tsv -fw -o negatives_pred_onehot.xlsx
