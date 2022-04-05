#!/usr/bin/env zsh
SRC="`git root`/src/chemistry"
$SRC/smiles2cid.R < negatives.smiles > negatives.tsv
mlr -I --tsv --from negatives.tsv filter '$cid != 0'
$SRC/pubchem_props.R -Hp smiles < negatives.tsv > negatives-pubchem_smiles.tsv
# check that the molecules found on pubchem are the same
cmp <($SRC/toCanonSmiles.py -H < negatives-pubchem_smiles.tsv) <(cut -f1 negatives.tsv | sed 1d) && rm negatives-pubchem_smiles.tsv
