#!/usr/bin/env zsh
# install muscle with conda install muscle -c bioconda
muscle -quiet < ../sequences.faa | fasta.py unwrap > muscle.faa
hmmbuild --amino muscle.{hmm,faa} > muscle.hmm.log
cat ../sequences.faa ../CAZy/seqs_len.faa | hmmalign --trim --amino --outformat A2M muscle.hmm - | fasta.py delete-lower > muscle.hmm.faa
fasta.py filter-quality -t 0.8 muscle.hmm.faa muscle_qual05.hmm.faa
