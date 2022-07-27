#!/usr/bin/env zsh
ROOT=`git root`
FEAT=`ls -d $ROOT/results/*-features/`
chemFeat=`ls $ROOT/results/*-chemicalFeatures/*_features.tsv`
blosum62Amb=$FEAT/blosum62Amb.tsv.gz
matchAmb=$FEAT/matchAmb.tsv.gz
train=$FEAT/train.tsv
testset=`ls $ROOT/results/*-experimentalValidation/experimentYield.tsv`

# randomforest.py -i enzyme cid -c reaction -n 10000 -o matchAmb.rf.joblib -a ${=chemFeat} $matchAmb < $train
# randomforest.py -i enzyme cid -c reaction -n 10000 -o blosum62Amb.rf.joblib -a ${=chemFeat} $blosum62Amb < $train

randomforest.py -m matchAmb.rf.joblib -a ${=chemFeat} $matchAmb < $testset > matchAmb-pred.tsv
randomforest.py -m blosum62Amb.rf.joblib -a ${=chemFeat} $blosum62Amb < $testset > blosum62Amb-pred.tsv

