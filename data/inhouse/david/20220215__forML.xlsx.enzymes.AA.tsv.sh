#!/usr/bin/env zsh
INFILE="20220215__forML.xlsx.enzymes.tsv"
OUTFILE="$INFILE:r.AA.tsv"

# get translation field of a file in EMBL format.
function EMBL_get_trans() {
    grep '^FT' | sed 's/FT *//' | tr -d '\n' | sed 's/.*translation="//' | sed 's/".*//'
}

echo "AA" > $OUTFILE

mlr --tsv --from $INFILE cut -f ENA | sed 1d |
while read ACC; do
    curl https://www.ebi.ac.uk/ena/browser/api/embl/$ACC | EMBL_get_trans >> $OUTFILE
    echo "" >> $OUTFILE # newline
done 

paste $INFILE $OUTFILE > temp && mv temp $OUTFILE
