#!/usr/bin/env zsh
# Generate a tsv with the train and test set for validating new experimental results.
# This would be training data with previously known reaction bool 
# and test set that is experimentally tested enzymes and acceptors vs each other and vs previous enzymes/acceptors.
ROOT=`git root`

TRAIN=`ls $ROOT/results/*features/rdkit-desc_muscle-ntern-hmm.tsv`

CID="$ROOT/data/inhouse/david/20220215__forML.xlsx.cas.tsv"
ENZ="$ROOT/data/inhouse/david/20220215__forML.xlsx.enzymes.AA.tsv"
acceptors=`ls $ROOT/results/*chemicalFeatures/acceptor_features.tsv`
alignment=`ls $ROOT/results/*align/muscle_qual05.hmm.nterm.tsv`

# subtable of enzyme features for the new ones we want to test
SEQS="seqs.tmp.tsv"
mlr --tsv --from $ENZ cut -f ENA then rename ENA,enzyme then\
    join -j enzyme -f $alignment > $SEQS

# subtable of acceptor features for the new ones we want to test
CHEMS="chems.tmp.tsv"
mlr --tsv --from $CID cut -f cid then join -j cid -f $acceptors > $CHEMS


ids="enzyme acceptor source cid reaction rate"

# for aaEnc in $ROOT/data/NCBI/{matchAmb,blosum62Amb}.tsv; do
for aaEnc in $ROOT/data/NCBI/matchAmb.tsv; do
    
    # encode_features.py -i ${=ids} --aa seq --aa-encoding "$aaEnc" \
    # "$acceptors" "$alignment" < $TRAIN |
    #     gzip > traintest_${aaEnc:r:t}.tsv.gz
    
    encode_features.py -i ${=ids} --aa seq --aa-encoding "$aaEnc" \
    "$CHEMS" "$alignment" < $TRAIN |
        gzip > traintest_chem_${aaEnc:r:t}.tsv.gz
    
    # encode_features.py -i ${=ids} --aa seq --aa-encoding "$aaEnc" \
    # "$acceptors" "$SEQS" < $TRAIN |
    #     gzip > traintest_seqs_${aaEnc:r:t}.tsv.gz
    
done

