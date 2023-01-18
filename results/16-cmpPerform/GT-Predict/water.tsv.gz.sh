#!/usr/bin/env zsh
# Run all-vs-all Smith-Waterman alignments with settings identical to default matlab swalign.
# EBLOSUM50 is identical to the blosum50 used by default by matlab, see
# https://au.mathworks.com/help/bioinfo/ref/blosum.html
# and the emboss/[share/EMBOSS/]data/EBLOSUM50 file.
ROOT=`git root`
GTPred=$ROOT/data/GT-Predict/Active_enzymes_protein_sequences.txt
SWopts='-sprotein1 -sprotein2 -datafile EBLOSUM50 -gapopen 8 -gapextend 8'
$ROOT/src/water.sh $GTPred ${=SWopts} < dataset1_unalign.fa | gzip > water.tsv.gz
