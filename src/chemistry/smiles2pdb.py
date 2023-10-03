#!/usr/bin/env python3
import argparse
import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolfiles

def get_parser():
    parser = argparse.ArgumentParser(description="Make .pdb files of structures from SMILES given in stdin.")
    parser.add_argument("id", nargs='?', default="cid", help="Column with filenames for the output PDBs. Default=\"cid\".")
    parser.add_argument("smiles", nargs='?', default="smiles", help="Name of column with SMILES strings.")
    parser.add_argument("-o", "--out", help="Output directory. Default=current working directory.")
    return parser


# example use
debug = False
if debug:
    args = get_parser().parse_args()
    df = pd.DataFrame(dict(cid=[40024], smiles=["CC1=CC2C(C(C1=O)O)(C3(CC(C(C34CO4)O2)O)C)CO"]))
else:
    args = get_parser().parse_args()
    df = pd.read_table(sys.stdin, sep='\t')

fnames = df[args.id]
smiles = df[args.smiles]
outdir = args.out

mols = [Chem.MolFromSmiles(s) for s in smiles]
mols = [Chem.AddHs(m) for m in mols]

# calculate 3D shapes
for i, m in enumerate(mols):
    AllChem.EmbedMolecule(m)
    while m.GetNumConformers() == 0:
        sys.stderr.write(f"Retrying CID={df.loc[i, args.id]} SMILES={df.loc[i, args.smiles]}\n")
        AllChem.EmbedMolecule(m)
sys.stderr.write(f"Conformers embedded.\n")


# make valid filenames
fnames = [str(f) for f in fnames]
if outdir: # will also work for outdir == ""
    fnames = [outdir + "/" + f for f in fnames]

for m, fname in zip(mols, fnames):
    if not fname.endswith('.pdb'): fname += '.pdb'
    Chem.rdmolfiles.MolToPDBFile(m, fname)
