#!/usr/bin/env zsh
[ -f ../pubchem/CID-Title.gz ] || ../pubchem/CID-Title.gz.sh

mlr -t --hi label cid,title ../pubchem/CID-Title.gz |
    mlr -t join -j cid -f rawAcceptor2cid.tsv > rawAcceptor_cid_title.tsv

mlr -t --from rawAcceptor2cid.tsv filter '$smiles != ""' + put '$title = $raw' > smiles.tmp.tsv

mlr -t template -f cid,smiles,raw,title rawAcceptor_cid_title.tsv smiles.tmp.tsv > temp && mv temp rawAcceptor_cid_title.tsv

nBefore=`cat rawAcceptor2cid.tsv | wc -l | xargs`
nAfter=`cat rawAcceptor_cid_title.tsv | wc -l | xargs`

# check that all titles are found
if [[ $nBefore != $nAfter ]]; then
    echo "$nBefore != $nAfter"
    exit 1
else
    rm smiles.tmp.tsv
fi

