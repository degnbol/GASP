#!/usr/bin/env zsh
# Seqs for 132 and 149 were provided by David. This checks that they are the same as what I already have.
# they are the same
mlr --tsv --from ForChristian__Results_screening_of_47_polyphenols_against_40_GTs.xlsx.enzymes.tsv filter '$ID == 132 || $ID == 149' then cut -f Accession,ID,AA |
    cmp - ForChristian__Results_screening_of_47_polyphenols_against_40_GTs.xlsx.enzymes-check_132+149.tsv
