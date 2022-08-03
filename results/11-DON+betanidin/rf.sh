#!/usr/bin/env zsh
ROOT=`git root`
chemFeat=`ls $ROOT/results/*-chemicalFeatures/{acceptors2,DON+betanidin}_features.tsv`
# sequence features
FEAT=`ls -d $ROOT/results/*-features`
blosum62Amb=$FEAT/blosum62Amb.tsv.gz
matchAmb=$FEAT/matchAmb.tsv.gz

# train
# randomforest.py -i enzyme cid -c reaction -n 10000 -o blosum62Amb.rf.joblib.gz -a ${=chemFeat} $blosum62Amb < DON+betanidin-trainset.tsv

# test
randomforest.py -m blosum62Amb.rf.joblib.gz -a ${=chemFeat} $blosum62Amb < DON-testset.tsv > DON-pred.tsv
randomforest.py -m blosum62Amb.rf.joblib.gz -a ${=chemFeat} $blosum62Amb < betanin-testset.tsv |
    mlr -t filter '$pred != ""' + remove-empty-columns > betaninCazy-pred.tsv
randomforest.py -m blosum62Amb.rf.joblib.gz -a ${=chemFeat} betanin_blosum62Amb.tsv.gz < betanin-testset.tsv |
    mlr -t filter '$pred != ""' > betanin-pred.tsv


