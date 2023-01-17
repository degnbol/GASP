#!/usr/bin/env zsh
ROOT=`git root`

INFILE=`\ls $ROOT/results/*generateNegatives/reactions.tsv`
SEQS=`\ls $ROOT/results/*-align/muscle_qual.hmm.nterm.tsv.gz`
CHEMS=`\ls $ROOT/results/*-chemicalFeatures/acceptors2_features.tsv`

# Trainset is GT-Predict data + negatives, and testset is dataset 1. So first we filter.
# annotate with unaligned sequences in "Enzymes" and smiles in "Metabolites".
# For simplicity, only keep unique entries for ESP required columns and mapping columns "enzyme", "cid".
# split in files of max 1000 entries as per limits of online predictor.
mlr --t2c --from $INFILE filter '$reaction != 0.5 && $source !=~ "^GT-Predict" && $source !=~ "^negatives"' +\
    join -j enzyme -f $SEQS + join -j cid --lk smiles -f $CHEMS +\
    put '$Enzymes = gsub($seq, "-", "")' +\
    rename smiles,Metabolites +\
    cut -f Enzymes,Metabolites,enzyme,cid + uniq -a +\
    split -n 1000 --prefix=ESP_input

echo "https://esp.cs.hhu.de/ES_pred_multiple"

