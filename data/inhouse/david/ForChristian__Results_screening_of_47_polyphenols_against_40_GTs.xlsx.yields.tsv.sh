#!/usr/bin/env zsh
# contents were copy-pasted then corrected for \r and % removed, then checked with mlr cat that there is no issues with windows chars.
FILE="ForChristian__Results_screening_of_47_polyphenols_against_40_GTs.xlsx.yields.tsv"
tr -d '%\r' < $FILE |
    mlr --tsv cat > temp && mv temp $FILE
