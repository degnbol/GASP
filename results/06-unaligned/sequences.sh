#!/usr/bin/env zsh
# Collect unaligned GT1 sequences.
ROOT=`git root`
DATA=$ROOT/data
LIB=$ROOT/tools/degnlib/subfunctions

# Initial GT-Predict, in-house
sed '/>/!s/\*//' "$DATA/GT-Predict/Active_enzymes_protein_sequences.txt" \
    "$DATA/inhouse/hts_ugt/hts_ugt_constructs.fa" |
    sed '/^>/s/MT_71G1/Mt_71G1/' | sed 's/>UGT/>At_/' | $LIB/fasta_table.py |
    cat <(mlr --tsv --from "$DATA/inhouse/pTMH.tsv" uniq -f enzyme,sequence) - |
    mlr --tsv uniq -a > sequences.tsv

# GT-Predict ext
mlr -t uniq -f Protein,AA "$DATA/GT-Predict/extensions.AA.tsv" | sed 1d >> sequences.tsv

# lit
mlr -t cut -f enzyme,sequence "$DATA/lit/sequences.tsv" | sed 1d >> sequences.tsv
mlr -t uniq -f Protein,AA "$DATA/lit/lit.tsv" | sed 1d >> sequences.tsv

# add new seqs that were tested experimentally
mlr -t --from $DATA/inhouse/david/20220215__forML.xlsx.enzymes.AA.tsv cut -f ENA,AA | sed 1d >> sequences.tsv

mlr -t -I uniq -a sequences.tsv
# muscle in the alignment process step only takes the first word, so we can't have spaces.
tr ' ' '-' < sequences.tsv > temp && mv temp sequences.tsv

# the following code returns nothing when there are no duplicate entries, which is best.
# Duplicates just means a bit of wasted effort in alignment, potentially.
echo "Duplicate entries:"
# use -D , so space isn't assumed to be a in-field sep
table_3D.R sequence -D '|' < sequences.tsv | cut -f1 | grep '|'

# also make fasta version
$LIB/table_fasta.py -d "\t" -H sequences.tsv > sequences.faa
