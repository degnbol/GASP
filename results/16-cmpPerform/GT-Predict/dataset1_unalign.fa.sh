#!/usr/bin/env zsh
ROOT=`git root`
TOOL=$ROOT/tools/degnlib/subfunctions
dataset1=`\ls $ROOT/results/*-dataset1/dataset1.tsv`
mlr -t --from $dataset1 cut -f enzyme,seq + put '$seq = gsub($seq, "-", "")' |
    $TOOL/table_fasta.py -H > dataset1_unalign.fa
