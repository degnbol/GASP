#!/usr/bin/env zsh
ROOT=`git root`
mlr --tsv uniq -f cid $ROOT/data/inhouse/david/20220215__forML.xlsx.cas.tsv > 20220215.cid
