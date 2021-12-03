#!/usr/bin/env zsh
# Get the SMILES from CIDs for reactive chemicals.
mlr --tsv --from `git root`/data/reactions.tsv filter '$reaction == 1' then uniq -f cid |
    `git root`/src/chemistry/pubchem_props.R -Hp smiles > positives.tsv

