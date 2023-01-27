#!/usr/bin/env zsh
ROOT=`git root`
FEAT=`ls -d $ROOT/results/*-features/`
chemFeat=`ls $ROOT/results/*-chemicalFeatures/acceptors2_features.tsv`
blosum62Amb=$FEAT/blosum62Amb.tsv.gz
# trained on all data available from dataset 1
model=`ls $ROOT/results/*-experimentalValidation/smpl_blosum62Amb.rf.joblib.gz`

randomforest.py -m $model -a $chemFeat $blosum62Amb < testset.tsv > blosum62Amb-pred.tsv
# sort
mlr -I -t --from blosum62Amb-pred.tsv sort -nr pred

