#!/usr/bin/env zsh
newSeqs="`git root`/data/inhouse/david/20220215__forML.xlsx.enzymes.tsv"
mlr --from $newSeqs rename ENA,enzyme + \
    join -f experimentYield.tsv -j enzyme + \
    uniq -f cid,enzyme,Yield > seqTest.tsv
