#!/usr/bin/env zsh
ROOT=`git root`
WORK=`ls -d $ROOT/results/*features`
reactions=`ls $ROOT/results/*generateNegatives/reactions.tsv`
mlr -t cut -x -f rate,acceptor,source + filter '$reaction != 0.5' $reactions > "$WORK/train.tsv"
