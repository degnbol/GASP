#!/usr/bin/env zsh
# Generate a lot of features for chemicals given their pubchem id (CID) or SMILES.
# USAGE: pipeline.sh [-t <integer>] [-m <integer>] INFILE.cid OUTFILE.tsv [PREVIOUS]
# - -t: is to set number of processes to run for parallelization, 
# e.g. set to 1 to disable parallelization.
# Default = all available minus 1 detected by calling nproc. If the call fails fallback to 1.
# - -m: Number of MDS features. Default = 12.
# - INFILE contains a CID or SMILES on each line.
# - OUTFILE.tsv is the generated features.
# - PREVIOUS is an optional name for a previous job, i.e. the name used as INFILE earlier.
# It is important to supply this in case of extending a dataset,
# so the MDS projected points are in the same space.
# A folder with intermediate files are also generated named INFILE-props
usage() { echo "Usage: $0 [-t <integer>] [-m <integer>] INFILE OUTFILE.tsv [PREVIOUS]" 1>&2; exit 1; }

while getopts ":s:p:" o; do
    case "${o}" in
        t)
            NPROC=${OPTARG}
            # must be integer
            [[ $NPROC =~ '^[0-9]+$' ]] || usage
            ;;
        m)
            MDS=${OPTARG}
            # must be integer
            [[ $MDS =~ '^[0-9]+$' ]] || usage
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

INFILE="$1"
OUTFILE="$2"
PREVIOUS="$3"

# mandatory arguments
if [ -z "$INFILE" ] || [ -z "$OUTFILE" ]; then
    usage
fi
# set defaults for optional arguments
if [ -z "$NPROC" ]; then
    NPROC=`nproc 2> /dev/null | tr -d '\n'`
    if [ -z "$NPROC" ]; then
        NPROC=1
    else
        # leave one processing unit for regular functioning
        NPROC=$((NPROC-1))
    fi
fi
echo "# Using $NPROC processing unit(s)"
if [ -z "$MDS" ]; then
    MDS=12
fi

if [ -n "$PREVIOUS" ]; then
    # set to the path of the fingerprints generated in the previous run
    PREVIOUS="${PREVIOUS:r}-props/E3FP.fpz" # :r = remove extension
fi

WORK="${INFILE:r}-props" # :r = remove extension
mkdir -p "$WORK"
{
	echo '#!/usr/bin/env zsh'
	echo '# This folder contains chemistry features that were generated with the call'
	echo "$0 $@"
} > "$WORK/README.sh"

# Other scripts are located relative to this file
SRC="$0:h"

# make sure we can exit the whole pipeline with a simple ctrl+c.
trap 'exit 1' SIGINT SIGTERM EXIT

# get SMILES and other pubchem listed properties that are always listed.
echo '# get SMILES from CIDs'
if [ -s "$WORK/pubchem.tsv" ]; then
    echo "# SKIP. Already exists: $WORK/pubchem.tsv"
else
    $SRC/pubchem_props.R < "$INFILE" > "$WORK/pubchem.tsv"
fi

echo '# get RDKit descriptors from SMILES'
if [ -s "$WORK/RDKitDescriptors.tsv" ]; then
    echo "SKIP. Already exists: $WORK/RDKitDescriptors.tsv"
else
    $SRC/rdkit-descriptors.py --nproc $NPROC --pdb "$WORK/PDBs" < "$WORK/pubchem.tsv" > "$WORK/RDKitDescriptors.tsv"
fi

# activate env to make sure pymol can be found
which pymol > /dev/null || {
    echo 'pymol not found in PATH. rerun after `conda activate GT`'
    exit 1
}

echo '# calculate areas from PDBs using PyMol'
if [ -s "$WORK/areas.tsv" ]; then
    echo "# SKIP. Already exists: $WORK/areas.tsv"
else
    N=`ls $WORK/PDBs/*.pdb | wc -l | xargs`
    if [ "$NPROC" -gt 1 ]; then
        $SRC/pymol_areas.sh $WORK/PDBs/*.pdb > "$WORK/areas.tsv" &
    else
        $SRC/pymol_areas.sh $WORK/PDBs/*.pdb > "$WORK/areas.tsv"
    fi
fi

echo '# make volumes from PDBs'
if [ -s "$WORK/volumes.tsv" ]; then
    echo "# SKIP. Already exists: $WORK/volumes.tsv"
else
    # NOTE that this only worked when searching for the processes by name, 
    # which can match any other process of the same name.
    trap 'pkill -f ProteinVolume_1.3; exit 1' SIGINT SIGTERM EXIT
    if [ "$NPROC" -gt 1 ]; then
        $SRC/ProteinVolume.sh "$WORK/PDBs" > "$WORK/volumes.tsv" &
    else
        $SRC/ProteinVolume.sh "$WORK/PDBs" > "$WORK/volumes.tsv"
    fi
fi

echo '# make E3FP fingerprints using SMILES'
# Don't check for skipping existing file since we read the outfile and see if it lacks entries.
N=`sed 1d $WORK/pubchem.tsv | wc -l | xargs`
# suppressing stderr with verbose warnings from RDKit conformer gen
# $SRC/e3fp-fprints.py -ct "$NPROC" "$WORK/E3FP.fpz" < "$WORK/pubchem.tsv" # 2> /dev/null | $SRC/../progress.sh $N

if [ -s "$WORK/E3FP_MDS.tsv" ]; then
    echo "# SKIP. Already exists: $WORK/E3FP_MDS.tsv"
else
    echo '# use the chemical fingerprints to generate MDS features'
    julia --startup-file=no --project=$SRC/../../ $SRC/E3FP_features.jl -ck $MDS $PREVIOUS "$WORK/E3FP.fpz" | mlr -t rename 'id,cid' > "$WORK/E3FP_MDS.tsv"
fi

if [ "$NPROC" -gt 1 ]; then
    # make sure all parallel bg processes finished
    wait
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

