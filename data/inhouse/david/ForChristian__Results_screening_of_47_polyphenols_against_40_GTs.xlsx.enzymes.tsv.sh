#!/usr/bin/env zsh
# enzyme table region from first sheet copy pasted then corrected and columns named.
FILE='ForChristian__Results_screening_of_47_polyphenols_against_40_GTs.xlsx.enzymes.tsv'

tr -d '\r' $FILE < | sed -E 's/^UGT ?//' | cat <(echo "ID\tDescription\tUGT\tAccession") - > temp && mv temp $FILE

# UGT==- means no UGT entry
mlr -I --tsv --from $FILE put 'if ($UGT == "-") {$UGT = ""}'

# there is an entry missing for UGT number 9
echo "9\t\t\tC0HFA0" >> $FILE

# change one EMBL id to uniprot accession
sed 's/ATO74555.1/A0A2D1N4Z8/' $FILE > temp && mv temp $FILE

# get seqs
echo "Accession\tAA" > enzymes.tmp.tsv
mlr --tsv cut -f Accession $FILE | sed 1d |
while read acc; do
    echo -n "$acc\t"
    curl https://www.uniprot.org/uniprot/$acc.fasta | sed 1d | tr -d '\n'
    echo ""
done >> enzymes.tmp.tsv

mlr --tsv join -j Accession -f $FILE enzymes.tmp.tsv > temp && mv temp $FILE
rm enzymes.tmp.tsv

