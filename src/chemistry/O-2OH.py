#!/usr/bin/env python3
import sys
import pandas as pd
from src.chemistry.chem_utils.rdkit_utils import *
from src.utilities.io_utils import git_root

# USAGE: src/chemistry/O-2OH.py < infile > outfile where infile has columns cid,smiles and outfile gets columns cid,smiles_OH


debug = False
if debug:
    df = pd.read_table(git_root("acceptors/acceptor_smiles.tsv"))
else:
    df = pd.read_table(sys.stdin)


mols = smiles2molecules(*(df[k] for k in ["smiles", "cid"]))
mols = mols[has_any_sub(mols, "C[O-]", "c[O-]")]
mols = replace_sub(mols, "C[O-]", "CO") # aliphatic
mols = replace_sub(mols, "c[O-]", "cO") # aromatic/phenolate

out = pd.DataFrame(dict(cid=get_names(mols), smiles_OH=get_smiles(mols)))
out.to_csv(sys.stdout, sep='\t', index=False)
