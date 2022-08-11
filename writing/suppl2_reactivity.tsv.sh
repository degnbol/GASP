#!/usr/bin/env zsh
ROOT=`git root`
reactions=`ls $ROOT/results/*-generateNegatives/reactions.tsv`
mlr -t --from $reactions \
    put '$rawType = $rate == "" ? "bool" : "rate"' +\
    put 'if ($source == "pTMH" || $source == "HTS_UGT") {$source = "Dataset 1"}' +\
    put 'if ($source == "GT-Predict extensions") {$source = "GT-Predict"}' +\
    put 'if ($source == "\"\"" || $source == "sebastian" || $source =~ "^http.*") {$source = "literature"}' +\
    uniq -c -f rawType,source,reaction + sort -f source,rawType -n reaction \
    > suppl2_reactivity.tsv
    
