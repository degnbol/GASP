#!/usr/bin/env zsh

for file in muscle{,.hmm,_qual05.hmm}.faa; do
    length=$(fasta.py length $file | sort -u)
    let length=length/2
    fasta.py range $file -r 1-$length | fasta.py table -i enzyme -s seq > ${file:r}.nterm.tsv
done

