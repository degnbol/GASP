#!/usr/bin/env zsh
# Write pymol_areas.tsv given filenames of PDB files
# REQUIREMENT: pymol installed and accessible with command `pymol`
# USAGE: pymol_areas.sh PDBs/*.pdb > pymol_areas.tsv
echo "cid\tSESA"
for fname in $@; do
	echo -n "${fname:t:r}\t"
	pymol $fname -Q -cmd 'get_area' | sed 's/.*: //' | sed 's/ .*//'
done
