#!/usr/bin/env zsh
`git root`/src/chemistry/pipeline.sh acceptors2{.cid,_features.tsv}
# you could also update the old files to match the MDS space
mlr -t --from acceptors2_features.tsv join -j cid -f <(echo cid | cat - acceptors.cid) > acceptors_features.tsv
mlr -t --from acceptors2_features.tsv join -j cid -f <(echo cid | cat - 20220215.cid) > 20220215_features.tsv
