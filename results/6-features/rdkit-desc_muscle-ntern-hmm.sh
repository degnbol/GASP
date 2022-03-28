#!/usr/bin/env zsh
ROOT=`git root`
WORK=`ls -d $ROOT/results/*features`
chemicals=`ls -d $ROOT/results/*chemicalFeatures/acceptor_features.tsv`
reactions=`ls -d $ROOT/results/*generateNegatives/reactions.tsv`
alignment="$ROOT/data/align/muscle_qual05.hmm.nterm.tsv"

# create a feature set for classification where we discard uncertain data (reaction == 0.5)
mlr --tsv --from "$chemicals" join -f "$reactions" -j cid then join -f "$alignment" -j enzyme then \
    cut -x -f smiles then filter '$reaction != 0.5' > "$WORK/rdkit-desc_muscle-ntern-hmm.tsv"

mlr --tsv filter '$source == "GT-Predict"' "$WORK/rdkit-desc_muscle-ntern-hmm.tsv" > "$WORK/rdkit-desc_muscle-ntern-hmm_gtpred.tsv"
mlr --tsv filter '$source != "GT-Predict"' "$WORK/rdkit-desc_muscle-ntern-hmm.tsv" > "$WORK/rdkit-desc_muscle-ntern-hmm_notgtpred.tsv"

