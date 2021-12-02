#!/usr/bin/env zsh
# use defaults from matlab
# declaring potein, probably making use of BLOSUM50
water_scores.sh -gapopen 8 -gapextend 8 -sprotein1 -sprotein2 raw/Active_enzymes_protein_sequences.txt < raw/hts_ugt_constructs.fa | mlr --tsv rename seq1,enzyme,seq2,gtpred_enzyme > water.tsv
