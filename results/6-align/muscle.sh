#!/usr/bin/env zsh
# install muscle with conda install muscle -c bioconda

DATA=`git root`/data
LIB=`git root`/tools/degnlib/subfunctions
INFILE=$DATA/unaligned/sequences.faa
CAZY=$DATA/CAZy/seqs_len.faa

# align main sequences with muscle. takes time.
muscle -align $INFILE -output muscle.faa.tmp
$LIB/fasta_unwrap.py muscle.faa.tmp > muscle.faa && rm muscle.faa.tmp

# write Hidden Markov Model from the alignment.
hmmbuild --amino muscle.{hmm,faa} > muscle.hmm.log

# align all sequences to the HMM, then remove non-consensus.
cat $INFILE $CAZY | hmmalign --trim --amino --outformat A2M muscle.hmm - |
    $LIB/fasta_delete_lower.py > muscle.hmm.faa

# remove low quality alignments
$LIB/fasta_filter_quality.py -t 0.8 muscle.hmm.faa muscle_qual05.hmm.faa
