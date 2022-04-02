#!/usr/bin/env zsh
LIB=`git root`/tools/degnlib/subfunctions

for file in muscle{,.hmm,_qual05.hmm}.faa; do
    length=$($LIB/fasta_length.py $file | sort -u)
    let length=length/2
    $LIB/fasta_range.py $file -r 1-$length | $LIB/fasta_table.py -i enzyme -s seq > ${file:r}.nterm.tsv
done

