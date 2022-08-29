#!/usr/bin/env zsh
ROOT=`git root`
WORK=`ls -d $ROOT/results/*features`
reactions=`ls $ROOT/results/*generateNegatives/reactions.tsv`
mlr -t filter '$reaction != 0.5' + cut -f enzyme,reaction,cid,source $reactions > "$WORK/train.tsv"
mlr -t filter '$source !=~ "negatives"' + cut -f enzyme,reaction,cid train.tsv > "$WORK/train_negativesSample.tsv"
mlr -t filter '$source =~ "negatives"' + sample -k 2000 + cut -f enzyme,reaction,cid train.tsv | sed 1d >> "$WORK/train_negativesSample.tsv"
