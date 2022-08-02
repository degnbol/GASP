#!/usr/bin/env zsh
ROOT=`git root`
FEAT=`ls -d $ROOT/results/*-features/`
chemFeat=`ls $ROOT/results/*-chemicalFeatures/acceptors2_features.tsv`
blosum62Amb=$FEAT/blosum62Amb.tsv.gz
matchAmb=$FEAT/matchAmb.tsv.gz
trainset=DON+betanidin-trainset.tsv
testset=DON+betanidin-testset.tsv

randomforest.py -i enzyme cid -c reaction -n 10000 -o blosum62Amb.rf.joblib.gz -a ${=chemFeat} $blosum62Amb < $train
randomforest.py -m blosum62Amb.rf.joblib.gz -a ${=chemFeat} $blosum62Amb < $testset > blosum62Amb-pred.tsv
auc -c Yield -t 50 -p pred blosum62Amb-pred.tsv -n

