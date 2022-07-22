#!/usr/bin/env zsh
mlr --tsv --from CAZy_GT1_structure.tsv cut -f Uniprot | sed 1d |
    while read uniprot; do
        wget https://rest.uniprot.org/uniprotkb/$uniprot.fasta
    done
