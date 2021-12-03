#!/usr/bin/env zsh
# grab the collected acceptor features relevant to the CIDs for the select few acceptors that we will focus on.
mlr --tsv --from acceptor_features.tsv join -j cid -f "`git root`/data/acceptorsOfInterest/acceptors-of-interest.tsv" > acceptorsOfInterest.tsv
