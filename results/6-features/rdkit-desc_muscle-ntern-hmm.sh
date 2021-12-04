#!/usr/bin/env zsh
ROOT=`git root`
OUT=$ROOT/results/*features
chemicals="$ROOT/results/*chemicalFeatures/acceptor_features.tsv"
reactions="$ROOT/results/*generateNegatives/reactions.tsv"
alignment="$ROOT/data/align/muscle_qual05.hmm.nterm.tsv"

# create a feature set for classification where we discard uncertain data (reaction == 0.5)
mlr --tsv --from "$chemicals" join -f "$reactions" -j cid then join -f "$alignment" -j enzyme then \
    cut -x -f smiles then filter '$reaction != 0.5' > "$OUT/rdkit-desc_muscle-ntern-hmm.tsv"

mlr --tsv filter '$source == "GT-Predict"' rdkit-desc_muscle-ntern-hmm.tsv > "$OUT/rdkit-desc_muscle-ntern-hmm_gtpred.tsv"
mlr --tsv filter '$source != "GT-Predict"' rdkit-desc_muscle-ntern-hmm.tsv > "$OUT/rdkit-desc_muscle-ntern-hmm_notgtpred.tsv"

