#!/usr/bin/env zsh
sed '/>/!s/\*//' "raw/gtpred/Active_enzymes_protein_sequences.txt" "raw/hts_ugt/hts_ugt_constructs.fa" | sed '/^>/s/MT_71G1/Mt_71G1/' | sed 's/>UGT/>At_/' |
    fasta.py table | cat <(mlr --tsv --from pTMH.tsv uniq -f enzyme,sequence) - | mlr --tsv uniq -a > sequences.tsv
mlr --tsv cut -f enzyme,sequence "lit/sebastian/sequences.tsv" | sed 1d >> sequences.tsv
# the following code returns nothing when there are no duplicate entries, which is best.
table_3D.R sequence < sequences.tsv | cut -f1 | grep ' '
table.py fasta -H sequences.tsv > sequences.faa

