#!/usr/bin/env zsh
ROOT=`git root`
ids="enzyme acceptor source cid"

echo 'cid' | cat - "$ROOT/data/inhouse/david/negatives.cid" |
    mlr --tsv join -j cid -f $ROOT/results/*chemicalFeatures/acceptor_features.tsv > negatives_props.tsv

mlr --tsv uniq -f enzyme,seq $ROOT/gt/features/rdkit-desc_muscle-ntern-hmm.tsv > train_seqs.tsv

for aaEnc in $ROOT/data/ncbi/{matchAmb,blosum62Amb}.tsv; do
    # train set:
    cat $ROOT/gt/features/rdkit-desc_muscle-ntern-hmm.tsv |
        encode_features.py -i ${=ids} reaction rate --aa seq --aa-encoding $aaEnc negatives_props.tsv train_seqs.tsv > test_negatives_${aaEnc:r:t}.tsv
done

rm negatives_props.tsv

