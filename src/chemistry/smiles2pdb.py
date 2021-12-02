#!/usr/bin/env python3
import argparse, sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolfiles

def get_parser():
    parser = argparse.ArgumentParser(description="Make .pdb files of structures from SMILES given in stdin.")
    parser.add_argument("id", help="Column with filenames for the output pdbs.")
    parser.add_argument("smiles", nargs='?', default="smiles", help="Name of column with SMILES strings.")
    return parser


# example use
debug = False
if debug:
    df = pd.read_table("~/biosustain/gt/rdkit/acceptor_smiles.tsv")
    fnames = df.cid
    smiles = df.smiles
else:
    args = get_parser().parse_args()
    df = pd.read_table(sys.stdin, sep='\t')
    fnames = df[args.id]
    smiles = df[args.smiles]


mols = [Chem.MolFromSmiles(s) for s in smiles]
mols = [Chem.AddHs(m) for m in mols]
# calculate 3D shapes
for m in mols: AllChem.EmbedMolecule(m)

for m, fname in zip(mols, fnames):
    fname = str(fname)
    if not fname.endswith('.pdb'): fname += '.pdb'
    Chem.rdmolfiles.MolToPDBFile(m, fname)
