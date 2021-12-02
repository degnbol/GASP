#!/usr/bin/env python3
import argparse, sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolfiles
from src.utilities.io_utils import git_root

def get_parser():
    parser = argparse.ArgumentParser(description="Make .pdb files of structures from SMILES given in stdin.")
    parser.add_argument("id", help="Column with filenames for the output PDBs.")
    parser.add_argument("smiles", nargs='?', default="smiles", help="Name of column with SMILES strings.")
    parser.add_argument("-o", "--out", help="Output directory. Default=current working directory.")
    return parser


# example use
debug = False
if debug:
    df = pd.read_table(git_root("data/pubchem.tsv"))
    fnames = df.cid
    smiles = df.smiles
    outdir = args.out
else:
    args = get_parser().parse_args()
    df = pd.read_table(sys.stdin, sep='\t')
    fnames = df[args.id]
    smiles = df[args.smiles]
    outdir = args.out


mols = [Chem.MolFromSmiles(s) for s in smiles]
mols = [Chem.AddHs(m) for m in mols]
# calculate 3D shapes
for m in mols: AllChem.EmbedMolecule(m)

# make valid filenames
fnames = [str(f) for f in fnames]
if outdir: # will also work for outdir == ""
    fnames = [outdir + "/" + f for f in fnames]

for m, fname in zip(mols, fnames):
    if not fname.endswith('.pdb'): fname += '.pdb'
    Chem.rdmolfiles.MolToPDBFile(m, fname)
