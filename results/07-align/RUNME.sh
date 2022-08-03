#!/usr/bin/env zsh
# INSTALL muscle with conda install muscle -c bioconda
ROOT=`git root`
LIB=$ROOT/tools/degnlib/subfunctions
export PYTHONPATH="$PYTHONPATH:$ROOT/tools/degnlib"

INFILE=`ls $ROOT/results/*unaligned/sequences.faa`
CAZY=$ROOT/data/CAZy/seqs_len.faa

# align main sequences with muscle. takes time.
muscle -align $INFILE -output muscle.faa.tmp
$LIB/fasta_unwrap.py muscle.faa.tmp > muscle.faa && rm muscle.faa.tmp

# write Hidden Markov Model from the alignment.
hmmbuild --amino muscle.{hmm,faa} > muscle.hmm.log

# align all sequences to the HMM
cat $INFILE $CAZY | ./hmmalign.sh
gzip -c muscle_qual.hmm.nterm.tsv > muscle_qual.hmm.nterm.tsv.gz

# check which training set sequences were filtered out
mlr --hi join --np --ul -j 1 -f <(grep '>' muscle.faa) <(grep '>' muscle_qual.hmm.faa) | sed 1d | cut -c2- > discarded.enz
if [ -s discarded.enz ]; then
    echo "The following enzymes were discarded due to gappy alignment (excluding bulk CAZy curated enzymes):"
    cat discarded.enz
fi

