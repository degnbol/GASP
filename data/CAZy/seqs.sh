#!/usr/bin/env zsh
`git root`/src/genbank2seqs.R < GT1.genbank > seqs.faa
fasta.py table seqs.faa -H | mlr --tsv filter 'strlen($seq) >= 300 && strlen($seq) <= 600' | table.py fasta -Hd $'\t' > seqs_len.faa
