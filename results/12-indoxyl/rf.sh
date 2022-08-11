#!/usr/bin/env zsh
ROOT=`git root`
FEAT=`ls -d $ROOT/results/*-features/`
chemFeat=`ls $ROOT/results/*-chemicalFeatures/acceptors2_features.tsv`
blosum62Amb=$FEAT/blosum62Amb.tsv.gz
testset=indoxyl-testset.tsv
# trained on all data available after adding the experimental validation data (aka. set 2)
model=`ls $ROOT/results/*-DON+betanidin/blosum62Amb.rf.joblib.gz`

randomforest.py -m $model -a $chemFeat $blosum62Amb < $testset > blosum62Amb-pred.tsv
# sort
mlr -I -t --from blosum62Amb-pred.tsv sort -nr pred

