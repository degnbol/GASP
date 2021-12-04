#!/usr/bin/env zsh
ROOT=`git root`

acceptorsOfInterest=$ROOT/results/*chemicalFeatures/acceptorsOfInterest.tsv
alignment="$ROOT/data/align/muscle_qual05.hmm.nterm.tsv"

# for aaEnc in $ROOT/data/NCBI/{matchAmb,blosum62Amb}.tsv; do
for aaEnc in $ROOT/data/NCBI/matchAmb.tsv; do
    cat rdkit-desc_muscle-ntern-hmm.tsv |
        encode_features.py -i enzyme acceptor source cid reaction rate --aa seq --aa-encoding "$aaEnc" "$acceptorsOfInterest" "$alignment" |
            gzip > traintest_${aaEnc:r:t}.tsv.gz
done

