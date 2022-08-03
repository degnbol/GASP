#!/usr/bin/env zsh
ROOT=`git root`
LIB=$ROOT/tools/degnlib/subfunctions

mlr -t uniq -f enzyme,AA $ROOT/data/inhouse/david/betaninProd.tsv |
    $LIB/table_fasta.py -H | $ROOT/results/*-align/hmmalign.sh

# -k to keep all features to make sure we have the ones present in the model.
cat muscle_qual.hmm.nterm.tsv |
    encode_features.py -k -i enzyme --aa seq --aa-encoding $ROOT/data/NCBI/blosum62Amb.tsv |
    gzip -c > betanin_blosum62Amb.tsv.gz

