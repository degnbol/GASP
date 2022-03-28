#!/usr/bin/env zsh
# Get the SMILES from CIDs for all chemicals that may be reactive or not.
DATA=`git root`/data
mlr --tsv uniq -f cid "$DATA/reactions/reactions.tsv" "$DATA/acceptorsOfInterest/acceptors-of-interest.tsv" | sed 1d |
    cat - "$DATA/inhouse/david/negatives.cid" | sort -u |
    `git root`/src/chemistry/pubchem_props.R -p smiles > acceptors.tsv

