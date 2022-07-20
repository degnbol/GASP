#!/usr/bin/env zsh
# Writes a copy of INFILE.tsv to INFILE.aa.tsv with amino acid sequences curated from EMBL,
# where INFILE.tsv must have EMBL accession ids in a column given as second arg, or by default "Accession".
# USE: acc2aa.sh INFILE.tsv [Accession column name]
INFILE=$1
ACC_COL=$2
[ -n "$ACC_COL" ] || ACC_COL="Accession"
OUTFILE="$INFILE:r.AA.tsv"

# get translation field of a file in EMBL format.
function EMBL_get_trans() {
    grep '^FT' | sed 's/FT *//' | tr -d '\n' | sed 's/.*translation="//' | sed 's/".*//'
}

echo "AA" > $OUTFILE

mlr --tsv --from $INFILE cut -f $ACC_COL | sed 1d |
while read ACC; do
    curl https://www.ebi.ac.uk/ena/browser/api/embl/$ACC | EMBL_get_trans >> $OUTFILE
    echo "" >> $OUTFILE # newline
done 

TMP=$RANDOM.tmp
paste $INFILE $OUTFILE > $TMP && mv $TMP $OUTFILE
