#!/usr/bin/env zsh
# Add negatives to the reaction data with the reaction bool = false.
mlr -t --from ../*validateAcceptors/reactions.tsv uniq -f enzyme | sed 1d |
while read enzyme; do
    mlr -t --from negatives_cid.tsv cut -f cid,replacement + put '$enzyme = "'$enzyme'"'
done | mlr --tsv filter '$cid != "cid"' + \
put '$reaction = 0; $rate = ""; $source = "negatives " . $replacement; $acceptor = ""' + \
cut -x -f replacement > reactions_neg.tsv

mlr -t unsparsify ../*validateAcceptors/reactions.tsv reactions_neg.tsv > reactions.tsv
rm reactions_neg.tsv

