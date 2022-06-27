#!/usr/bin/env zsh
LIB=`git root`/tools/degnlib/subfunctions
$LIB/fasta_table.py -i enzyme -s seq gtpred.muscle.hmm.fa > gtpred.muscle.hmm.tsv
