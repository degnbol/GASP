#!/usr/bin/env zsh
# join experimental data with sequence data to get a
# yield table with at least columns cid,enzyme,yield

ROOT=`git root`
LIB="$ROOT/tools/degnlib/subfunctions"

expBasename="$ROOT/data/inhouse/david/ForChristian__Results_screening_of_47_polyphenols_against_40_GTs.xlsx"
expEnz="$expBasename.enzymes.tsv"
expYield="$expBasename.tsv"
trainEnz=`ls $ROOT/results/*unaligned/sequences.faa`
CAZY="$ROOT/data/CAZy/seqs_len.faa"

# match on sequence.
# Take first match from $trainEnz $CAZY cat which means we have preference for $trainEnz matches.
cat $trainEnz $CAZY | $LIB/fasta_table.py -i enzyme -s AA |
    mlr --tsv join -j AA -f $expEnz then reorder -e -f AA then\
    sort -n ID then head -n 1 -g ID > experimentEnz.tsv
# 43 lines but expEnz has 47 lines. These 4 are unmatched:
cat $trainEnz $CAZY | $LIB/fasta_table.py -i enzyme -s AA |
    mlr --tsv join -j AA -f $expEnz --np --ul > experimentEnz-unmatched.tsv

# can join on Accession, or ID, it was checked that they produce identical tables (after mlr reorder).
mlr --tsv --from experimentEnz.tsv cut -x -f AA then\
    join -j Accession -f $expYield > experimentYield.tsv




