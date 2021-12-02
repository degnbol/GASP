#!/usr/bin/env zsh
echo $'cid\ttitle' | cat - ../../data/ncbi/pubchem/raw/CID-Title | mlr --tsvlite join -j cid -f raw-acceptor2cid.tsv > raw-acceptor_cid_title.tsv
