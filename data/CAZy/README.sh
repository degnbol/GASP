#!/usr/bin/env zsh
# cazy parser was modified to fix a bug and and only use GT1
./create_cazy_db.py.sh
# then it was run, which crawls the cazy website
./create_cazy_db.py
# then a small grep -v step which results in the table cazy.tsv with header: name, species, genbank.
./cazy.tsv.sh
# unique genbank ids are collected
./GT1.genbank.sh
# use rentrez service to find amino acid sequences for all the genbank IDs (seqs.faa),
# then discard sequences shorter than 300 AAs or longer than 600 (seqs_len.faa).
./seqs.sh
