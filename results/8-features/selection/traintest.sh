#!/usr/bin/env zsh
ROOT=`git root`

# for mat in match blosum62; do
for mat in match; do
    cat $ROOT/gt/features/rdkit-desc_muscle-ntern-hmm.tsv |
        mlr --tsv put 'if (substr($enzyme, 0, 1) == "At") {$set = "At"} elif (substr($enzyme, 0, 1) == "Zm") {$set = "Zm"} else {$set = "others"}' |
    encode_features.py -i enzyme acceptor source cid reaction rate set --aa seq --aa-encoding $ROOT/data/ncbi/$mat.tsv |
    gzip > traintest_$mat.tsv.gz
done
