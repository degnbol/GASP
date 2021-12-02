#!/usr/bin/env python3
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys
from rdkit.Chem.AtomPairs import Pairs


df = pd.read_table("/gt/acceptors/acceptor_smiles.tsv")
smiles = df.smiles

mols = [Chem.MolFromSmiles(s) for s in smiles]

# is calculating 3D shapes necessary?
mols = [Chem.AddHs(m) for m in mols]
for m in mols: AllChem.EmbedMolecule(m)

# fingerprints listed here:
# http://www.rdkit.org/docs/GettingStartedInPython.html#list-of-available-fingerprints
rdkfps = [Chem.RDKFingerprint(m) for m in mols]
maccs = [MACCSkeys.GenMACCSKeys(m) for m in mols]
pairfps = [Pairs.GetAtomPairFingerprint(m) for m in mols]
morgan = [AllChem.GetMorganFingerprint(m, 2) for m in mols] # second arg is int radius, using 2 since it was used in e3fp manual

# TODO use one or more of the features to get pairwise distances between molecules and from that the points in abstract feature space



