#!/usr/bin/env zsh
ROOT=`git root`
LIB=$ROOT/tools/degnlib/subfunctions
$LIB/table_excel.py blosum62Amb-pred.tsv -s 'indoxyl vs all seqs' -o blosum62Amb-pred.xlsx -fw
