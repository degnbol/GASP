#!/usr/bin/env zsh
ROOT=`git root`
FEAT=`ls -d $ROOT/results/*-features/`
# chemFeat=`ls $ROOT/results/*-chemicalFeatures/{20220215,acceptors}_features.tsv`
chemFeat=`ls $ROOT/results/*-chemicalFeatures/acceptors2_features.tsv`
blosum62Amb=$FEAT/blosum62Amb.tsv.gz
matchAmb=$FEAT/matchAmb.tsv.gz
train=$FEAT/train.tsv
trainSmpl=$FEAT/train_negativesSample.tsv
testsetAll=experimentYield.tsv
testsetDist=testsetDist.tsv
testsetSeq=testsetNovelSeqs.tsv
testsetChem=testsetNovelChems.tsv

# randomforest.py -i enzyme cid -c reaction -n 1000 -o matchAmb.rf.joblib.gz -a ${=chemFeat} $matchAmb < $train
# randomforest.py -i enzyme cid -c reaction -n 1000 -o blosum62Amb.rf.joblib.gz -a ${=chemFeat} $blosum62Amb < $train
# randomforest.py -i enzyme cid -c reaction -n 1000 --class-weight balanced -o balance_matchAmb.rf.joblib.gz -a ${=chemFeat} $matchAmb < $train
# randomforest.py -i enzyme cid -c reaction -n 1000 --class-weight balanced -o balance_blosum62Amb.rf.joblib.gz -a ${=chemFeat} $blosum62Amb < $train
# randomforest.py -i enzyme cid -c reaction -n 1000 -o smpl_matchAmb.rf.joblib.gz -a ${=chemFeat} $matchAmb < $trainSmpl
# randomforest.py -i enzyme cid -c reaction -n 1000 -o smpl_blosum62Amb.rf.joblib.gz -a ${=chemFeat} $blosum62Amb < $trainSmpl
# randomforest.py -i enzyme cid -c reaction -n 1000 --class-weight balanced -o smpl_balance_matchAmb.rf.joblib.gz -a ${=chemFeat} $matchAmb < $trainSmpl
# randomforest.py -i enzyme cid -c reaction -n 1000 --class-weight balanced -o smpl_balance_blosum62Amb.rf.joblib.gz -a ${=chemFeat} $blosum62Amb < $trainSmpl

# echo match
# for model in *matchAmb.rf.joblib.gz; do
#     echo $model
#     echo all
#     randomforest.py -m $model -a ${=chemFeat} $matchAmb < $testsetAll  | auc -c Yield -t 50 -p pred -n -a
#     echo dist
#     randomforest.py -m $model -a ${=chemFeat} $matchAmb < $testsetDist | auc -c Yield -t 50 -p pred -n -a
#     echo seq
#     randomforest.py -m $model -a ${=chemFeat} $matchAmb < $testsetSeq  | auc -c Yield -t 50 -p pred -n -a
#     echo chem
#     randomforest.py -m $model -a ${=chemFeat} $matchAmb < $testsetChem | auc -c Yield -t 50 -p pred -n -a
# done
#
# echo blosum
# for model in *blosum62Amb.rf.joblib.gz; do
#     echo $model
#     echo all
#     randomforest.py -m $model -a ${=chemFeat} $blosum62Amb < $testsetAll  | auc -c Yield -t 50 -p pred -n -a
#     echo dist
#     randomforest.py -m $model -a ${=chemFeat} $blosum62Amb < $testsetDist | auc -c Yield -t 50 -p pred -n -a
#     echo seq
#     randomforest.py -m $model -a ${=chemFeat} $blosum62Amb < $testsetSeq  | auc -c Yield -t 50 -p pred -n -a
#     echo chem
#     randomforest.py -m $model -a ${=chemFeat} $blosum62Amb < $testsetChem | auc -c Yield -t 50 -p pred -n -a
# done

# randomforest.py -m smpl_matchAmb.rf.joblib.gz -a ${=chemFeat} $matchAmb < $testsetAll > smpl_matchAmb-pred.tsv
# auc smpl_matchAmb-pred.tsv -c Yield -t 50 -p pred -n -a

randomforest.py -m smpl_blosum62Amb.rf.joblib.gz -a ${=chemFeat} $blosum62Amb < $testsetAll > smpl_blosum62Amb-pred.tsv
auc smpl_blosum62Amb-pred.tsv -c Yield -t 50 -p pred -n -a
