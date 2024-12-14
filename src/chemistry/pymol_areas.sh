#!/usr/bin/env zsh
# Write pymol_areas.tsv given filenames of PDB files
# REQUIREMENT: pymol installed and accessible with command `pymol`
# USAGE: pymol_areas.sh PDBs/*.pdb > pymol_areas.tsv

# -k = don't load ~/.pymolrc which may affect things and cause logging outputs
# -c = command-line only mode
# -q = suppress startup messages
# -Q = suppress all texts
pymol -kcxqQ ${=@} $0:h/get_area.pml.py

