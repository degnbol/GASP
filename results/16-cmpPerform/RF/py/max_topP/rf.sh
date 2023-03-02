#!/usr/bin/env zsh
ROOT=`git root`
chemFeat=`\ls $ROOT/results/*-chemicalFeatures/acceptors2_features.tsv`
blosum62Amb=`\ls -d $ROOT/results/*-features/blosum62Amb.tsv.gz`

featNames=`tr '\n' ' ' < $ROOT/results/*-negFeatSelect/feature_selection-topP.colname`

# train
randomforest.py -i enzyme cid -c reaction -n 1000 -o blosum62Amb.rf.joblib.gz -a $chemFeat $blosum62Amb -F ${=featNames} -R 'seq_*' < ../trainset.tsv
# test
randomforest.py -i enzyme cid -c reaction -m blosum62Amb.rf.joblib.gz -a ${=chemFeat} $blosum62Amb -F ${=featNames} -R 'seq_*' < ../testset.tsv > pred_blosum62Amb.tsv
# AUC
mlr -t filter '$pred != ""' pred_blosum62Amb.tsv | auc -c reaction -p pred > auc.tsv
rocauc.py -c reaction -p pred pred_blosum62Amb.tsv -o roc-topP.png --title topP
