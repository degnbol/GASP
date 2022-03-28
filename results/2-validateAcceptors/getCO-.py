#!/usr/bin/env python3
# Find acceptors with O- (base) that should be converted to OH (acid). 
# Can also be used to assert that non are found in the validated set of acceptors.
# USAGE: ./getCO-.py < INFILE.tsv > OUTFILE
# INFILE has columns smiles,cid. OUTFILE has a single column named cid
from src.chemistry.chem_utils.rdkit_utils import *
import pandas as pd
import sys

# cid,smiles
acceptors = pd.read_table(sys.stdin)

mols = smiles2molecules(*(acceptors[k] for k in ["smiles", "cid"]))
bases = mols[has_any_sub(mols, "C[O-]", "c[O-]")]

print("cid")
for base in bases:
    print(Molecule(base).name)

