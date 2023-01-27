#!/usr/bin/env zsh
ROOT=`git root`

INFILE=`\ls $ROOT/results/*generateNegatives/reactions.tsv`
SEQS=`\ls $ROOT/results/*-align/muscle_qual.hmm.nterm.tsv.gz`
CHEMS=`\ls $ROOT/results/*-chemicalFeatures/acceptors2_features.tsv`

# Trainset is GT-Predict data + negatives (k=2000), and testset is dataset 1. So first we filter.
mlr -t --from $INFILE filter '$reaction != 0.5 && $source =~ "^GT-Predict"' +\
    join -j enzyme -f $SEQS + join -j cid --lk smiles -f $CHEMS + uniq -a > trainset.tsv
mlr -t --from $INFILE filter '$reaction != 0.5 && $source =~ "^negatives"'  +\
    join -j enzyme -f $SEQS + join -j cid --lk smiles -f $CHEMS + uniq -a +\
    sample -k 2000 | sed 1d >> trainset.tsv

# rate is always missing
mlr -t -I cut -x -f rate trainset.tsv
