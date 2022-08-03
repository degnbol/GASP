#!/usr/bin/env zsh
# Align sequences to hmm built from training data sequences in results/*-align/,
# silently filter out low quality (>25% gap) alignments, then take N-term half.
# USE: results/*-align/hmmalign.sh < INFILE.faa
# WRITES: muscle.hmm.a2m, muscle.hmm.faa, and muscle_qual.hmm.faa, muscle_qual.hmm.nterm.tsv in pwd
ROOT=`git root`
LIB=$ROOT/tools/degnlib/subfunctions
export PYTHONPATH="$PYTHONPATH:$ROOT/tools/degnlib"

# align all sequences to the HMM
hmmalign --trim --amino --outformat A2M $ROOT/results/*-align/muscle.hmm - > muscle.hmm.a2m
# remove non-consensus
$LIB/fasta_delete_lower.py muscle.hmm.a2m > muscle.hmm.faa
# remove low quality alignments
$LIB/fasta_filter_quality.py -t 0.75 muscle.hmm.faa muscle_qual.hmm.faa

# use only N-term (discard last half of alignments)
length=$($LIB/fasta_length.py muscle_qual.hmm.faa | sort -u)
let length=length/2
$LIB/fasta_range.py muscle_qual.hmm.faa -r 1-$length |
    $LIB/fasta_table.py -i enzyme -s seq > muscle_qual.hmm.nterm.tsv

