#!/usr/bin/env zsh
# header
echo "Protein\tAccession\tAA" > extensions.AA.tsv

mlr --tsv --from CAZy_GT1_structure.tsv cut -f 'Protein Name,Uniprot' | sed 1d |
    while read protein uniprot; do
        echo -n "$protein\t$uniprot\t"
        sed 1d $uniprot.fasta | tr -d '\n'
        echo '' # newline
    done >> extensions.AA.tsv

LIB=`git root`/tools/degnlib/subfunctions
cat *.faa | $LIB/fasta_table.py | tr ' ' '\t' | cut -f2- >> extensions.AA.tsv
