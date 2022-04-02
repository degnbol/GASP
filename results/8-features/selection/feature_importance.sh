#!/usr/bin/env zsh
cd $HOME/biosustain/gt/feature_importance
ids="enzyme acceptor cid source"
chemFeatures=$(zcat traintest_match.tsv.gz | head -n1 | tr '\t' '\n' | grep -v seq_ | tr '\n' ' ')
zcat traintest_match.tsv.gz | ./feature_importance.py -c reaction -i ${=ids} -S ${=chemFeatures} -s set -e rate -t 1 --nproc 30 > feature_importance.log
