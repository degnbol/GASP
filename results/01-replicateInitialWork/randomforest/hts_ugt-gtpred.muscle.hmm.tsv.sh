#!/usr/bin/env zsh
LIB="`git root`/tools/degnlib/subfunctions"
$LIB/fasta_table.py -i enzyme -s seq hts_ugt-gtpred.muscle.hmm.fa > hts_ugt-gtpred.muscle.hmm.tsv
