#!/usr/bin/env zsh
export PYTHONPATH="$PYTHONPATH:`git root`/tools/degnlib"
LIB=`git root`/tools/degnlib/subfunctions

cat `git root`/data/GT-Predict/Active_enzymes_protein_sequences.txt | muscle -quiet 2> gtpred.muscle.log |
    $LIB/fasta_unwrap.py > gtpred.muscle.fa
