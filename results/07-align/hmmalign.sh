#!/usr/bin/env zsh
# Align sequences to hmm built from training data sequences in results/*-align/.
# Low quality (>25% gap) alignments will be silently filtered out.
# USE: results/*-align/hmmalign.sh < INFILE.faa
# WRITES: muscle.hmm.a2m, muscle.hmm.faa, and muscle_qual.hmm.faa in pwd
ROOT=`git root`
LIB=$ROOT/tools/degnlib/subfunctions
export PYTHONPATH="$PYTHONPATH:$ROOT/tools/degnlib"

# align all sequences to the HMM
hmmalign --trim --amino --outformat A2M $ROOT/results/*-align/muscle.hmm - > muscle.hmm.a2m
# remove non-consensus
$LIB/fasta_delete_lower.py muscle.hmm.a2m > muscle.hmm.faa
# remove low quality alignments
$LIB/fasta_filter_quality.py -t 0.75 muscle.hmm.faa muscle_qual.hmm.faa

