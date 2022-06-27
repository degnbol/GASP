#!/usr/bin/env zsh
LIB=`git root`/tools/degnlib/subfunctions
hmmalign --trim --amino --outformat A2M gtpred.muscle.hmm `git root`/data/GT-Predict/Active_enzymes_protein_sequences.txt |
    $LIB/fasta_delete_lower.py > gtpred.muscle.hmm.fa
