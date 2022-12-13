#!/usr/bin/env zsh
ROOT=`git root`
LIB=$ROOT/tools/degnlib/subfunctions

# mlr -t --from blosum62Amb-pred.tsv join -j enzyme -f ../07-align/muscle_qual_wUnalign.hmm.nterm.tsv.gz > blosum62Amb-pred-wSeqs.tsv
mlr -t --from ../*-experimentalValidation/experimentEnz.tsv join -j enzyme -f blosum62Amb-pred-wSeqs.tsv > blosum62Amb-pred-wSeqs-dataset2.tsv

$LIB/table_excel.py blosum62Amb-pred-wSeqs{,-dataset2}.tsv -s 'aloesone vs all seqs' 'aloesone vs dataset 2' -o blosum62Amb-pred.xlsx -fw

