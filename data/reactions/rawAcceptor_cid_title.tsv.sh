#!/usr/bin/env zsh
# you may need to first run `cd ../pubchem; ./CID-Title.sh; cd -`
echo $'cid\ttitle' | cat - ../pubchem/CID-Title |
    mlr --tsv join -j cid -f rawAcceptor2cid.tsv > rawAcceptor_cid_title.tsv
