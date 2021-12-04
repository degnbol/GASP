#!/usr/bin/env zsh
# Purpose: replicate results of GT-Predict. 
# They simply use nearest neighbor for the prediction given an enzyme. We use Smith-Waterman alignment to find the closest matching enzyme from GT-Predict measured enzymes given any query enzyme.
# The query enzymes that we test are enzymes in-house enzymes from the Enzyme Engineering group.
#
# use defaults from GT-Predict MATLAB code to exactly replicate their approach. 
# declaring potein, probably making use of BLOSUM50
water_scores.sh -gapopen 8 -gapextend 8 -sprotein1 -sprotein2 raw/Active_enzymes_protein_sequences.txt < raw/hts_ugt_constructs.fa | mlr --tsv rename seq1,enzyme,seq2,gtpred_enzyme > water.tsv
