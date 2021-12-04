#!/usr/bin/env zsh
# Add negatives to the reaction data with the reaction bool = false.
mlr --tsv --from `git root`/data/reactions.tsv uniq -f enzyme | sed 1d | while read enzyme; do
    mlr --tsv --from `git root`/acceptors/gen_negatives/negatives.tsv cut -f cid then \
        put '$enzyme = "'$enzyme'"'
done | mlr --tsv filter '$cid != "cid"' then \
    put '$reaction = 0; $rate = ""; $source = "negatives"; $acceptor = ""' |
    mlrcat.sh `git root`/data/reactions.tsv - > `git root`/data/reactions-neg.tsv
