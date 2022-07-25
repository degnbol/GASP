#!/usr/bin/env zsh
# Get the SMILES from CIDs for all chemicals that may be reactive or not.
DATA=`git root`/data
mlr --tsv uniq -f cid "$DATA/reactions/reactions.tsv" \
    "$DATA/acceptorsOfInterest/acceptors-of-interest.tsv" \
    "$DATA/inhouse/david/20220215__forML.xlsx.cas.tsv" | sed 1d | sort -u |
    `git root`/src/chemistry/pubchem_props.R -p smiles |
    ./getCO-.py > CO-.tsv

