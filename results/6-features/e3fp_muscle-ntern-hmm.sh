#!/usr/bin/env zsh
E3FP_MDS=../*chemicalFeatures/acceptors-props/E3FP_MDS.tsv
reactions=../*generateNegatives/reactions.tsv
alignment="`git root`/data/align/muscle_qual05.hmm.nterm.tsv"

mlr --tsv --from $E3FP_MDS join -f "$reactions" -j cid then join -f "$alignment" -j enzyme then cut -x -f smiles then filter '$reaction != 0.5' > e3fp_muscle-ntern-hmm.tsv

#mlr --tsv filter '$source == "GT-Predict"' e3fp_muscle-ntern-hmm.tsv > e3fp_muscle-ntern-hmm_gtpred.tsv
#mlr --tsv filter '$source != "GT-Predict"' e3fp_muscle-ntern-hmm.tsv > e3fp_muscle-ntern-hmm_notgtpred.tsv


