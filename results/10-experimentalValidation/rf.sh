#!/usr/bin/env zsh
ROOT=`git root`
FEAT=`ls -d $ROOT/results/*-features/`
# chemFeat=`ls $ROOT/results/*-chemicalFeatures/{20220215,acceptors}_features.tsv`
chemFeat=`ls $ROOT/results/*-chemicalFeatures/acceptors2_features.tsv`
blosum62Amb=$FEAT/blosum62Amb.tsv.gz
matchAmb=$FEAT/matchAmb.tsv.gz
train=$FEAT/train.tsv
testset=seqTest.tsv
testsetAll=`ls $ROOT/results/*-experimentalValidation/experimentYield.tsv`

randomforest.py -i enzyme cid -c reaction -n 10000 -o matchAmb.rf.joblib.gz -a ${=chemFeat} $matchAmb < $train
randomforest.py -i enzyme cid -c reaction -n 10000 -o blosum62Amb.rf.joblib.gz -a ${=chemFeat} $blosum62Amb < $train

randomforest.py -m matchAmb.rf.joblib.gz -a ${=chemFeat} $matchAmb < $testset > matchAmb-pred.tsv
randomforest.py -m blosum62Amb.rf.joblib.gz -a ${=chemFeat} $blosum62Amb < $testset > blosum62Amb-pred.tsv
randomforest.py -m blosum62Amb.rf.joblib.gz -a ${=chemFeat} $blosum62Amb < $testsetAll > pred_all.tsv

auc -c Yield -t 50 -p pred matchAmb-pred.tsv -n
auc -c Yield -t 50 -p pred blosum62Amb-pred.tsv -n
auc -c Yield -t 50 -p pred pred_all.tsv -n

