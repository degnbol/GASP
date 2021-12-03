#!/usr/bin/env zsh
# use conda env gt that has rdkit
# To add new acceptors, edit acceptors.cid.sh and categories/categories.sh
echo '# collect all cids'
./acceptors.cid.sh
echo '# generate negatives from positives'
./gen_negatives/commands.sh

`git root`/src/chemistry/pipeline.sh acceptors.cid acceptor_features.tsv

echo '# extract entries for the acceptors of interest'
cd acceptors-of-interest
./acceptor_props.tsv.sh
cd ..
