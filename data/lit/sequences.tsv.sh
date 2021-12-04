#!/usr/bin/env zsh
mlr --icsv --otsv --from `git root`/data/lit/raw/sebastian/gt_all_data.csv rename 'C-GT,enzyme,aa_sequence,sequence,Uniprot/Genbank,uniprot/genbank' then \
    uniq -f enzyme,sequence,isoelectric_point,molecular_weight,uniprot/genbank > sequences.tsv
