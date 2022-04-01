#!/usr/bin/env zsh
# Generate a lot of features of chemicals given their pubchem id (CID).
# USAGE: pipeline.sh INFILE OUTFILE.tsv
# A folder with intermediate files are also generated named INFILE-props
INFILE="$1"
OUTFILE="$2"

WORK="${INFILE:r}-props" # remove extension
mkdir -p "$WORK"
{
	echo '#!/usr/bin/env zsh'
	echo '# This folder contains chemistry features that were generated with the call'
	echo "$0 $@"
} > "$WORK/README.sh"

SRC="`git root`/src/chemistry"

# get SMILES and other pubchem listed properties that are always listed.
echo '# get SMILES from CIDs'
$SRC/pubchem_props.R < "$INFILE" > "$WORK/pubchem.tsv"

echo '# get RDKit descriptors from SMILES'
$SRC/rdkit-descriptors.py < "$WORK/pubchem.tsv" > "$WORK/RDKitDescriptors.tsv"

echo '# make PDBs from SMILES'
mkdir -p "$WORK/PDBs"
$SRC/smiles2pdb.py cid -o "$WORK/PDBs" < "$WORK/pubchem.tsv"

echo '# calculate areas from PDBs using PyMol'
N=`ls $WORK/PDBs/*.pdb | wc -l | xargs`
$SRC/pymol_areas.sh $WORK/PDBs/*.pdb "$WORK/areas.tsv" | $SRC/../progress.sh $N

echo '# make volumes from PDBs'
$SRC/ProteinVolume.sh "$WORK/PDBs" > "$WORK/volumes.tsv"

echo '# make E3FP fingerprints using SMILES'
$SRC/e3fp-fprints.py "$WORK/E3FP.fpz" < "$WORK/pubchem.tsv"

echo '# use the compressed fingerprints to generate 12 MDS features'
$SRC/e3fp-features.py "$WORK/E3FP.fpz" -ck 12 | mlr --tsv rename 'id,cid' > "$WORK/E3FP_MDS.tsv"

echo '# join intermediate files'
mlr --tsv   --from "$WORK/volumes.tsv" \
    join -j cid -f "$WORK/RDKitDescriptors.tsv" then \
    join -j cid -f "$WORK/areas.tsv" then \
    join -j cid -f "$WORK/E3FP_MDS.tsv" then \
    rename -r ' ,_' > "$OUTFILE"

