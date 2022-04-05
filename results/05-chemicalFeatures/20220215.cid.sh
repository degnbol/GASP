#!/usr/bin/env zsh
CID="`git root`/data/inhouse/david/20220215__forML.xlsx.cas.tsv"
mlr --tsv uniq -f cid $CID | sed 1d > 20220215.cid
