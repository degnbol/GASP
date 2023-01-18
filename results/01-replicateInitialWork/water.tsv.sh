#!/usr/bin/env zsh
# Purpose: replicate results of GT-Predict. 
# They simply use nearest neighbor for the prediction given an enzyme.
# They use Smith-Waterman alignment to find the closest matching enzyme from GT-Predict measured enzymes given any query enzyme.
# They simply call matlab swalign(seq1, seq2) which by default is:
# - for amino acids, which by default uses BLOSUM50.
# - gap open penalty = 8
# - gap extend penalty = gap open penalty
# https://au.mathworks.com/help/bioinfo/ref/swalign.html

# The query enzymes that we test are enzymes in-house enzymes from the Enzyme Engineering group.

water_scores.sh -gapopen 8 -gapextend 8 -sprotein1 -sprotein2 raw/Active_enzymes_protein_sequences.txt < raw/hts_ugt_constructs.fa | mlr --tsv rename seq1,enzyme,seq2,gtpred_enzyme > water.tsv
