#!/usr/bin/env zsh
# not used but should give the same result.
ROOT=`git root`
testset=`\ls $ROOT/results/*-dataset1/dataset1.tsv`
CHEM=`\ls $ROOT/results/*-chemicalFeatures/acceptors2_features.tsv`
$ROOT/src/classifier.jl trainset.tsv $CHEM -t $testset -j cid enzyme -c reaction rate -o RF.jld -m n_estimators=1000 -s seq blosum62Amb 
