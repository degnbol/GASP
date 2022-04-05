#!/usr/bin/env zsh
# enzyme table region from first sheet copy pasted then corrected and columns named.
FILE='ForChristian__Results_screening_of_47_polyphenols_against_40_GTs.xlsx.enzymes.tsv'
tr -d '\r' $FILE < | sed -E 's/^UGT ?//' | cat <(echo "ID\tDescription\tUGT\tAccession") - > temp && mv temp $FILE
