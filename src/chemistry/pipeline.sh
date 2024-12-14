#!/usr/bin/env zsh
# Generate a lot of features for chemicals given their pubchem id (CID).
# USAGE: pipeline.sh INFILE.cid OUTFILE.tsv [PREVIOUS]
# INFILE contains a CID on each line.
# OUTFILE.tsv is the generated features.
# PREVIOUS is an optional name for a previous job, i.e. the name used as INFILE earlier.
# It is important to supply this if new chemicals are added to a previous run,
# so the MDS projected points are in the same space.
# A folder with intermediate files are also generated named INFILE-props
INFILE="$1"
OUTFILE="$2"
PREVIOUS="$3"

WORK="${INFILE:r}-props" # :r = remove extension
if [ -n "$PREVIOUS" ]; then
    # set to the path of the fingerprints generated in the previous run
    PREVIOUS="${PREVIOUS:r}-props/E3FP.fpz"
fi

mkdir -p "$WORK"
{
	echo '#!/usr/bin/env zsh'
	echo '# This folder contains chemistry features that were generated with the call'
	echo "$0 $@"
} > "$WORK/README.sh"

SRC="$0:h"

# get SMILES and other pubchem listed properties that are always listed.
echo '# get SMILES from CIDs'
if [ -s "$WORK/pubchem.tsv" ]; then
    echo "SKIP. Already exists: $WORK/pubchem.tsv"
else
    $SRC/pubchem_props.R < "$INFILE" > "$WORK/pubchem.tsv"
fi

echo '# get RDKit descriptors from SMILES'
if [ -s "$WORK/RDKitDescriptors.tsv" ]; then
    echo "SKIP. Already exists: $WORK/RDKitDescriptors.tsv"
else
    $SRC/rdkit-descriptors.py --pdb "$WORK/PDBs" < "$WORK/pubchem.tsv" > "$WORK/RDKitDescriptors.tsv"
fi

# activate env to make sure pymol can be found
which pymol > /dev/null || {
    echo 'pymol not found in PATH. rerun after `conda activate GT`'
    exit 1
}

echo '# calculate areas from PDBs using PyMol'
if [ -s "$WORK/areas.tsv" ]; then
    echo "SKIP. Already exists: $WORK/areas.tsv"
else
    N=`ls $WORK/PDBs/*.pdb | wc -l | xargs`
    $SRC/pymol_areas.sh $WORK/PDBs/*.pdb > "$WORK/areas.tsv"
fi

echo '# make volumes from PDBs'
if [ -s "$WORK/volumes.tsv" ]; then
    echo "SKIP. Already exists: $WORK/volumes.tsv"
else
    $SRC/ProteinVolume.sh "$WORK/PDBs" > "$WORK/volumes.tsv"
fi

echo '# make E3FP fingerprints using SMILES'
if [ -s "$WORK/E3FP.fpz" ]; then
    echo "SKIP. Already exists: $WORK/E3FP.fpz"
else
    # because of python subprocess weirdness we had to print to stderr,
    # then redirect pipe to stdout with 2>&1 in order to use progress.sh
    # grep for numbers since there can be other stderr warnings. Note that these are silenced by the progress reporting.
    $SRC/e3fp-fprints.py "$WORK/E3FP.fpz" -c < "$WORK/pubchem.tsv" 2>&1 | grep -E '^[0-9]+$' | $SRC/../progress.sh $N
fi

if [ -s "$WORK/E3FP_MDS.tsv" ]; then
    echo "SKIP. Already exists: $WORK/E3FP_MDS.tsv"
else
    echo '# use the chemical fingerprints to generate 12 MDS features'
    $SRC/E3FP_features.jl -ck 12 $PREVIOUS "$WORK/E3FP.fpz" | mlr -t rename 'id,cid' > "$WORK/E3FP_MDS.tsv"
fi

# RDKit and ProteinVolume both generate a Van Der Waals Volume.
# We keep the version from ProteinVolume to have consistency between different volume calculations.
# This is done by having it earlier in the list, since miller takes priority to earlier files in the join.
echo '# join intermediate files'
mlr -t      --from "$WORK/volumes.tsv" \
    join -j cid -f "$WORK/RDKitDescriptors.tsv" + \
    join -j cid -f "$WORK/areas.tsv" + \
    join -j cid -f "$WORK/E3FP_MDS.tsv" + \
    rename -r ' ,_' > "$OUTFILE"

