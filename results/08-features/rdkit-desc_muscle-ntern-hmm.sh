#!/usr/bin/env zsh
ROOT=`git root`
WORK=`ls -d $ROOT/results/*features`
chemicals=`ls $ROOT/results/*chemicalFeatures/acceptor_features.tsv`
reactions=`ls $ROOT/results/*generateNegatives/reactions.tsv`
alignment=`ls $ROOT/results/*align/muscle_qual.hmm.nterm.tsv`

# create a feature set for classification where we discard uncertain data (reaction == 0.5)
mlr --tsv --from "$chemicals" join -f "$reactions" -j cid then join -f "$alignment" -j enzyme then \
    cut -x -f smiles,rate then filter '$reaction != 0.5' > "$WORK/rdkit-desc_muscle-ntern-hmm.tsv"

gzip -c "$WORK/rdkit-desc_muscle-ntern-hmm.tsv" >  "$WORK/rdkit-desc_muscle-ntern-hmm.tsv.gz"

mlr --tsv filter '$source == "GT-Predict"' "$WORK/rdkit-desc_muscle-ntern-hmm.tsv" > "$WORK/rdkit-desc_muscle-ntern-hmm_gtpred.tsv"
mlr --tsv filter '$source != "GT-Predict"' "$WORK/rdkit-desc_muscle-ntern-hmm.tsv" > "$WORK/rdkit-desc_muscle-ntern-hmm_notgtpred.tsv"

