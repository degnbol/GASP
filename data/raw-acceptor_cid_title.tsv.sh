#!/usr/bin/env zsh
# you may need to first run `cd pubchem; ./CID-Title.sh; cd ..`
echo $'cid\ttitle' | cat - pubchem/CID-Title | mlr --tsvlite join -j cid -f raw-acceptor2cid.tsv > raw-acceptor_cid_title.tsv
