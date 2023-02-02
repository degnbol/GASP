#!/usr/bin/env zsh
ROOT=`git root`
chemFeat=`\ls $ROOT/results/*-chemicalFeatures/acceptors2_features.tsv`
blosum62Amb=`\ls -d $ROOT/results/*-features/blosum62Amb.tsv.gz`

reactions=`\ls $ROOT/results/*generateNegatives/reactions.tsv`
# train on GT-Predict and a sample of negatives, test on everythin else.
mlr -t filter '$reaction != 0.5 && $source =~ "GT-Predict"' + cut -x -f rate,acceptor,source $reactions > trainset.tsv
# mlr -t filter '$reaction != 0.5 && $source =~ "negatives"' + sample -k 2000 + cut -x -f rate,acceptor,source $reactions | sed 1d >> trainset.tsv
mlr -t filter '$reaction != 0.5 && ($source =~ "HTS_UGT" || $source =~ "pTMH")' + cut -x -f rate,acceptor,source $reactions > testset.tsv

# train
randomforest.py -i enzyme cid -c reaction -n 1000 -o blosum62Amb.rf.joblib.gz -a $chemFeat $blosum62Amb < trainset.tsv
# test
randomforest.py -m blosum62Amb.rf.joblib.gz -a ${=chemFeat} $blosum62Amb < testset.tsv > pred_blosum62Amb.tsv
# AUC
mlr -t filter '$pred != ""' pred_blosum62Amb.tsv | auc -c reaction -p pred > auc.tsv
rocauc.py -c reaction -p pred pred_blosum62Amb.tsv -o roc.png --title allFeats
