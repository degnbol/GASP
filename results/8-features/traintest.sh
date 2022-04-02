#!/usr/bin/env zsh
ROOT=`git root`
for aaEnc in $ROOT/data/NCBI/{match,blosum62}.tsv; do
cat rdkit-desc_muscle-ntern-hmm_gtpred.tsv |
    encode_features.py -i enzyme acceptor source cid reaction rate --aa seq --aa-encoding $aaEnc rdkit-desc_muscle-ntern-hmm_notgtpred.tsv |
        gzip > traintest_${aaEnc:r:t}.tsv.gz
done
