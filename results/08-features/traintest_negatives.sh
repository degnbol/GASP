#!/usr/bin/env zsh
ROOT=`git root`
ids="enzyme acceptor source cid"

mlr --tsv --from $ROOT/results/*-generateNegatives/negatives.tsv uniq -f cid then \
    join -j cid -f $ROOT/results/*chemicalFeatures/acceptor_features.tsv > negatives_props.tsv

mlr --tsv uniq -f enzyme,seq $ROOT/results/*features/rdkit-desc_muscle-ntern-hmm.tsv > train_seqs.tsv

for aaEnc in $ROOT/data/NCBI/{matchAmb,blosum62Amb}.tsv; do
    # train set:
    cat $ROOT/results/*features/rdkit-desc_muscle-ntern-hmm.tsv |
        encode_features.py -i ${=ids} reaction rate --aa seq --aa-encoding $aaEnc negatives_props.tsv train_seqs.tsv > test_negatives_${aaEnc:r:t}.tsv
done

rm negatives_props.tsv

