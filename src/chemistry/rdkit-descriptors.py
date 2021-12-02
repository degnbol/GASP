#!/usr/bin/env python3
import argparse, sys
import pandas as pd
import src.chemistry.chem_utils.rdkit_utils as rd
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, rdFreeSASA

# descriptors are listed here:
# https://www.rdkit.org/docs/GettingStartedInPython.html#list-of-available-descriptors

def get_parser():
    parser = argparse.ArgumentParser(
        usage="rdkit-descriptors.py < infile.tsv > outfile.tsv",
        description="Calculate physio-chemical properties from SMILES using RDKit.")
    parser.add_argument("smiles", nargs='?', default='smiles', help="Name of column with SMILES.")
    return parser


debug = False
if debug:
    df = pd.read_table("~/biosustain/gt/acceptors/acceptor_smiles.tsv")
    smiles = df.smiles
else:
    args = get_parser().parse_args()
    df = pd.read_table(sys.stdin)
    smiles = df[args.smiles]


mols = [Chem.MolFromSmiles(s) for s in smiles]
mols = [Chem.AddHs(m) for m in mols]
# calculate 3D shapes
for m in mols: AllChem.EmbedMolecule(m)

funcs_pubchem = ["MolWt", "NumHAcceptors", "NumHDonors", "NumRotatableBonds"]

funcs = ["NumValenceElectrons", "NumHeteroatoms",
         "MinPartialCharge", "MaxPartialCharge", "MinAbsPartialCharge", "MaxAbsPartialCharge"]
for f in funcs:
    df[f] = [getattr(Descriptors, f)(m) for m in mols]

funcs = ["CalcHallKierAlpha", "CalcLabuteASA",
         "CalcNumRings", "CalcNumAmideBonds", "CalcNumBridgeheadAtoms", "CalcNumHeterocycles",
         "CalcNumLipinskiHBA", "CalcNumLipinskiHBD", "CalcNumSpiroAtoms",
         # "CalcTPSA", # this is also found in PubChem
         "CalcNumAliphaticRings", "CalcNumAliphaticCarbocycles", "CalcNumAliphaticHeterocycles",
         "CalcNumAromaticRings", "CalcNumAromaticCarbocycles", "CalcNumAromaticHeterocycles",
         "CalcNumSaturatedCarbocycles", "CalcNumSaturatedHeterocycles", "CalcNumSaturatedRings",
         "CalcFractionCSP3"]
# requires a 3D conformer
funcs_3D = ["CalcEccentricity", "CalcAsphericity", "CalcInertialShapeFactor", "CalcNPR1", "CalcNPR2", "CalcSpherocityIndex",
            "CalcRadiusOfGyration", ]
# error from these:
funcs_err = ["CalcNumAtomStereoCenters"]
# descriptors that are not scalar. Currently not used.
funcs_array = ["BCUT2D"]

for f in funcs + funcs_3D:
    colname = f
    if colname.startswith('Calc'): colname = colname[len('Calc'):]
    df[colname] = [getattr(Descriptors.rdMolDescriptors, f)(m) for m in mols]

df["MolLogP"], df["MolMR"] = zip(*[Descriptors.rdMolDescriptors.CalcCrippenDescriptors(m) for m in mols])
df["HeavyAtomCount"] = [m.GetNumHeavyAtoms() for m in mols]
df["VanDerWaalsVolume"] = [AllChem.ComputeMolVolume(m) for m in mols]
df["SASA"] = [rdFreeSASA.CalcSASA(m, rdFreeSASA.classifyAtoms(m)) for m in mols]

# from http://rdkit.org/docs/source/rdkit.Chem.Fragments.html
for name, f in rd.numFragments.items(): df[name] = [f(m) for m in mols]

for name in ['primary_alcohol', 'secondary_alcohol', 'tertiary_alcohol', 'coumarin', 'thiophenol', 'hydroxyamino']:
    df[name] = [rd.get_num_substructs(m, rd.substructs[name]) for m in mols]


# df.to_csv("~/biosustain/gt/acceptors/rdkit-descriptors.tsv.tmp", sep='\t', index=False)
df.to_csv(sys.stdout, sep='\t', index=False)
