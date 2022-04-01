#!/usr/bin/env zsh
# Collect unaligned GT1 sequences.

DATA=`git root`/data
LIB=`git root`/tools/degnlib/subfunctions

sed '/>/!s/\*//' "$DATA/GT-Predict/Active_enzymes_protein_sequences.txt" \
    "$DATA/inhouse/hts_ugt/hts_ugt_constructs.fa" |
    sed '/^>/s/MT_71G1/Mt_71G1/' | sed 's/>UGT/>At_/' | $LIB/fasta_table.py |
    cat <(mlr --tsv --from "$DATA/inhouse/pTMH.tsv" uniq -f enzyme,sequence) - |
    mlr --tsv uniq -a > sequences.tsv
mlr --tsv cut -f enzyme,sequence "$DATA/lit/sequences.tsv" | sed 1d >> sequences.tsv

# add new seqs that were tested experimentally
mlr --tsv --from $DATA/inhouse/david/20220215__forML.xlsx.enzymes.AA.tsv cut -f ENA,AA | sed 1d >> sequences.tsv

# the following code returns nothing when there are no duplicate entries, which is best.
echo "Duplicate entries:"
table_3D.R sequence < sequences.tsv | cut -f1 | grep ' '

# also make a fasta copy
$LIB/table_fasta.py -H sequences.tsv > sequences.faa

