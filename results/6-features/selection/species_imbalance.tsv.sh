#!/usr/bin/env zsh
cut -c-2 $ROOT/gt/features/rdkit-desc_muscle-ntern-hmm.tsv | mlr --tsv uniq -c -f en > species_imbalance.tsv
