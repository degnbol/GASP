#!/usr/bin/env zsh
for file in pred*.tsv; do
    echo -n $file:r | sed $'s/-yield/\t/'
    mlr --tsv --from $file stats2 -a corr -f Yield,pred | sed 1d
done > corr.tsv
