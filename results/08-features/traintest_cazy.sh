#!/usr/bin/env zsh
ROOT=`git root`

acceptorsOfInterest=`ls $ROOT/results/*chemicalFeatures/acceptorsOfInterest.tsv`
alignment=`ls $ROOT/results/*align/muscle_qual05.hmm.nterm.tsv`

ids="enzyme acceptor source cid reaction rate"

# for aaEnc in $ROOT/data/NCBI/{matchAmb,blosum62Amb}.tsv; do
for aaEnc in $ROOT/data/NCBI/matchAmb.tsv; do
    cat rdkit-desc_muscle-ntern-hmm.tsv |
        encode_features.py -i ${=ids} --aa seq --aa-encoding "$aaEnc" \
        "$acceptorsOfInterest" "$alignment" |
            gzip > traintest_${aaEnc:r:t}.tsv.gz
done

