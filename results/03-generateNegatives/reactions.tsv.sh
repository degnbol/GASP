#!/usr/bin/env zsh
LIB=`git root`/tools/degnlib/

# Add negatives to the reaction data with the reaction bool = false.
mlr --tsv --from ../*validateAcceptors/reactions.tsv uniq -f enzyme | sed 1d |
while read enzyme; do
    mlr --tsv --from negatives.tsv cut -f cid then put '$enzyme = "'$enzyme'"'
done |
    mlr --tsv filter '$cid != "cid"' then\
    put '$reaction = 0; $rate = ""; $source = "negatives"; $acceptor = ""' |
    $LIB/mlr-cat ../*validateAcceptors/reactions.tsv /dev/stdin > reactions.tsv
