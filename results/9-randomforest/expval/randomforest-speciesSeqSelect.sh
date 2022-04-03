#!/usr/bin/env zsh
ROOT=`git root`
INFILE=traintest_chem_matchAmb.tsv.gz
chemFeatures=$(zcat $INFILE | head -n1 | tr '\t' '\n' | grep -v 'seq_' | tr '\n' ' ')
features=$(cat $ROOT/results/*features/selection/species_select.txt | tr '\n' ' ')
ids="enzyme acceptor source cid"
zcat $INFILE | randomforest.py -t 1 -i ${=ids} -c reaction -e rate -n 10000 -F ${=chemFeatures} ${=features} > pred_chem_matchAmb-speciesSeqSelect.tsv

