#!/usr/bin/env zsh
export PYTHONPATH="$PYTHONPATH:`git root`/tools/degnlib"
LIB=`git root`/tools/degnlib/subfunctions

infile="`git root`/data/GT-Predict/Active_enzymes_protein_sequences.txt"
muscle -align $infile -output gtpred.muscle.fa -quiet
$LIB/fasta_unwrap.py gtpred.muscle.fa temp && mv temp gtpred.muscle.fa
