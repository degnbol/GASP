#!/usr/bin/env zsh
# grab the collected acceptor features relevant to the CIDs for the select few acceptors that we will focus on.
aoi=`git root`/data/acceptorsOfInterest/acceptors-of-interest.tsv
mlr -t --from acceptors_features.tsv join -j cid -f $aoi > acceptorsOfInterest.tsv
