#!/usr/bin/env zsh
ROOT=`git root`
for file in $ROOT/results/*randomforest/expval/pred*.tsv.gz; do
    mlr --tsv --from $file join -j cid,enzyme -f experimentYield.tsv > ${file:t:r:r}-yield.tsv
done
