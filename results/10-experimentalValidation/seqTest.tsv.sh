#!/usr/bin/env zsh
newSeqs="`git root`/data/inhouse/david/20220215__forML.xlsx.enzymes.tsv"
# match on UGT_code instead instead ENA vs enzyme since we then would lose 2/8 
# enzymes due to no match between different kinds of namings/ids.
mlr --from $newSeqs rename UGT_code,ID + \
    join -f experimentYield.tsv -j ID + \
    uniq -f cid,enzyme,Yield > seqTest.tsv
