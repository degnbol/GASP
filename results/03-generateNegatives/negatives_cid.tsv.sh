#!/usr/bin/env zsh
SRC="`git root`/src/chemistry"
mlr -t --from negatives.tsv uniq -f smiles | sed 1d | $SRC/smiles2cid.R |
    mlr -t filter '$cid != 0' + join -j smiles -f negatives.tsv > negatives_cid.tsv

echo "checking that the molecules found on pubchem are the same..."
mlr -t uniq -f cid negatives_cid.tsv | sed 1d | $SRC/pubchem_props.R |
    $SRC/toCanonSmiles.py -Ho pubchemSmiles |
    mlr -t cut -f pubchemSmiles,cid |
    mlr -t join -j cid -f negatives_cid.tsv > negatives_cid_pubchemSmiles.tsv

mlr -t filter '$smiles != $pubchemSmiles' negatives_cid_pubchemSmiles.tsv > negatives_cid_pubchemSmiles-diff.tsv
echo "Printing any entries with differences:"
cat negatives_cid_pubchemSmiles.tsv

if [[ `wc -l negatives_cid_pubchemSmiles-diff.tsv` -eq 0 ]]; then
    rm negatives_cid_pubchemSmiles{,-diff}.tsv
fi

