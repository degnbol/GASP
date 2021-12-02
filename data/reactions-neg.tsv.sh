#!/usr/bin/env zsh
. setroot.zsh
mlr --tsv --from $ROOT/gt/data/reactions.tsv uniq -f enzyme | sed 1d | while read enzyme; do
mlr --tsv --from $ROOT/gt/acceptors/gen_negatives/negatives.tsv cut -f cid then put '$enzyme = "'$enzyme'"'; done |
mlr --tsv filter '$cid != "cid"' then put '$reaction = 0; $rate = ""; $source = "negatives"; $acceptor = ""' |
mlrcat.sh $ROOT/gt/data/reactions.tsv - > $ROOT/gt/data/reactions-neg.tsv
