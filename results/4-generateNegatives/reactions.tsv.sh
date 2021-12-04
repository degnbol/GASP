#!/usr/bin/env zsh
# Add negatives to the reaction data with the reaction bool = false.
mlr --tsv --from ../*validateAcceptors/reactions.tsv uniq -f enzyme | sed 1d | while read enzyme; do
    mlr --tsv --from negatives.tsv cut -f cid then put '$enzyme = "'$enzyme'"'
done | mlr --tsv filter '$cid != "cid"' then put '$reaction = 0; $rate = ""; $source = "negatives"; $acceptor = ""' |
    mlrcat.sh ../*validateAcceptors/reactions.tsv - > reactions.tsv
