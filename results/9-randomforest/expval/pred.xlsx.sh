#!/usr/bin/env zsh
ROOT=`git root`

# for file in pred*.gz; do
#     gzcat $file | mlr --tsv cut -f cid,enzyme,pred > $file:r
# done

grep '^seq' $ROOT/results/*features/selection/species_select.txt > speciesSeqSelect.tmp
INFILES="pred.xlsx.info speciesSeqSelect.tmp pred_chem_matchAmb-speciesSeqSelect.tsv pred_chem_matchAmb.tsv pred_seqs_matchAmb-speciesSeqSelect.tsv pred_seqs_matchAmb.tsv"
table.py excel ${=INFILES} -s "Overview" "species selected seq features" "chem match speciesSeqSelect" "chem match" "seqs match speciesSeqSelect" "seqs match" -o "pred.xlsx" -f


