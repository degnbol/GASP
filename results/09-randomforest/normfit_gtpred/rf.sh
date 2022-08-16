#!/usr/bin/env zsh
ROOT=`git root`
chemFeat=`ls $ROOT/results/*-chemicalFeatures/acceptors_features.tsv`
# sequence features
FEAT=`ls -d $ROOT/results/*-features`
matchAmb=$FEAT/matchAmb.tsv.gz
blosum62Amb=$FEAT/blosum62Amb.tsv.gz

reactions=`ls $ROOT/results/*generateNegatives/reactions.tsv`
# train on GT-Predict, test on everything else. Ignoring negatives for now.
mlr -t filter '$reaction != 0.5 && ($source =~ "GT-Predict" && $source != "negatives")' + cut -x -f rate,acceptor,source $reactions > trainset.tsv
mlr -t filter '$reaction != 0.5 && ($source !=~ "GT-Predict" && $source != "negatives")' + cut -x -f rate,acceptor,source $reactions > testset.tsv

# train
randomforest.py -i enzyme cid -c reaction -n 10000 -o matchAmb.rf.joblib.gz -a $chemFeat $matchAmb < trainset.tsv
randomforest.py -i enzyme cid -c reaction -n 10000 -o blosum62Amb.rf.joblib.gz -a $chemFeat $blosum62Amb < trainset.tsv

# test
randomforest.py -m matchAmb.rf.joblib.gz -a ${=chemFeat} $matchAmb < testset.tsv > pred_matchAmb.tsv
randomforest.py -m blosum62Amb.rf.joblib.gz -a ${=chemFeat} $blosum62Amb < testset.tsv > pred_blosum62Amb.tsv

# AUC
mlr -t filter '$pred != ""' pred_matchAmb.tsv | auc -c reaction -p pred
mlr -t filter '$pred != ""' pred_blosum62Amb.tsv | auc -c reaction -p pred
