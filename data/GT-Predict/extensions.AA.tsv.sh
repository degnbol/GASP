#!/usr/bin/env zsh
# remove species indication in protein name.
mlr --tsv --from extensions.tsv put '$Protein = sub($Protein, "^Mt", "")' then \
    join -j Protein -f seqs/extensions.AA.tsv > extensions.AA.tsv
