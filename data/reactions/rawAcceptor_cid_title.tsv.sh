#!/usr/bin/env zsh
[ -f ../pubchem/CID-Title.gz ] || ../pubchem/CID-Title.gz.sh

# mlr -t --hi label cid,title ../pubchem/CID-Title.gz |
#     mlr -t join -j cid -f rawAcceptor2cid.tsv > rawAcceptor_cid_title.tsv

nBefore=`cat rawAcceptor2cid.tsv | wc -l | xargs`
nAfter=`cat rawAcceptor_cid_title.tsv | wc -l | xargs`

# check that all titles are found
if [[ $nBefore != $nAfter ]]; then
    echo "$nBefore != $nAfter"
    exit 1
fi

