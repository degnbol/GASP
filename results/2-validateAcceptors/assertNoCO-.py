#!/usr/bin/env python3
from src.chemistry.chem_utils.rdkit_utils import *
from src.utilities.io_utils import git_root

# SMILES from chemicals
acceptors = pd.read_table("acceptors.tsv")

mols = smiles2molecules(*(acceptors[k] for k in ["smiles", "cid"]))
assert not has_any_sub(mols, "C[O-]", "c[O-]").any(), "There are some acceptors with O- that should be converted to OH."

