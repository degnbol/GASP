#!/usr/bin/env zsh
# Write pymol_areas.tsv given filenames of PDB files
# REQUIREMENT: pymol installed and accessible with command `pymol`
# USAGE: pymol_areas.sh PDBs/*.pdb pymol_areas.tsv > LOG
INFILES=${@:1:-1}
OUTFILE=$@[-1]

# make sure OUTFILE is supplied.
if [[ "${OUTFILE:e}" == "pdb" ]]; then
	echo "OUTFILE not supplied" 1>&2
	exit 1
fi

# after pymol 2.5 they started ignoring atoms with flag ignore,
# which apparently was on all atoms
CMD='flag ignore, all, clear; get_area' 

echo "cid\tSESA" > $OUTFILE
for fname in ${=INFILES}; do
	name=${fname:t:r}
	echo $name
	echo -n "$name\t" >> $OUTFILE
	pymol $fname -Q -cmd $CMD | grep -v '^PyMOL>' | grep -oE ' [0-9.]+' | tr -d ' ' >> $OUTFILE
done
