#!/usr/bin/env zsh
# install muscle with conda install muscle -c bioconda
export PYTHONPATH="$PYTHONPATH:`git root`/tools/degnlib"

ROOT=`git root`
LIB=$ROOT/tools/degnlib/subfunctions
INFILE=`ls $ROOT/results/*unaligned/sequences.faa`
CAZY=$ROOT/data/CAZy/seqs_len.faa

# align main sequences with muscle. takes time.
muscle -align $INFILE -output muscle.faa.tmp
$LIB/fasta_unwrap.py muscle.faa.tmp > muscle.faa && rm muscle.faa.tmp

# write Hidden Markov Model from the alignment.
hmmbuild --amino muscle.{hmm,faa} > muscle.hmm.log

# align all sequences to the HMM
cat $INFILE $CAZY | hmmalign --trim --amino --outformat A2M muscle.hmm - > muscle.hmm.a2m
# remove non-consensus
$LIB/fasta_delete_lower.py muscle.hmm.a2m > muscle.hmm.faa

# remove low quality alignments
$LIB/fasta_filter_quality.py -t 0.75 muscle.hmm.faa muscle_qual.hmm.faa

mlr --hi join --np --ul -j 1 -f <(grep '>' muscle.faa) <(grep '>' muscle_qual.hmm.faa) | sed 1d | cut -c2- > discarded.enz
if [ -s discarded.enz ]; then
    echo "The following enzymes were discarded due to gappy alignment (excluding bulk CAZy curated enzymes):"
    cat discarded.enz
fi

