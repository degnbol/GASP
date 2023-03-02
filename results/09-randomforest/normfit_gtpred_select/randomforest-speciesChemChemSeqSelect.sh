#!/usr/bin/env zsh
ROOT=`git root`
FEATURES=`\ls -d $ROOT/results/*features`
ids="enzyme acceptor source cid"
chemFeatures=$(cat $FEATURES/selection/species_select.txt | grep -v 'seq_' | tr '\n' ' ')
seqFeatures=$(cat $FEATURES/selection/selected_ratethres.txt | grep 'seq_' | tr '\n' ' ')
gzcat $FEATURES/traintest_match.tsv | randomforest.py -t 1 -i ${=ids} -c reaction -e rate -n 10000 -F ${=chemFeatures} ${=seqFeatures} > pred_onehot-speciesChemChemSeqSelect.tsv
gzcat $FEATURES/traintest_blosum62.tsv | randomforest.py -t 1 -i ${=ids} -c reaction -e rate -n 10000 -F ${=chemFeatures} ${=seqFeatures} > pred_blosum-speciesChemChemSeqSelect.tsv
