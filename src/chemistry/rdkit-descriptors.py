#!/usr/bin/env python3
import argparse
import sys, os
import numpy as np
import pandas as pd
import chem_utils.rdkit_utils as rd
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, rdFreeSASA
from multiprocess import Pool

# descriptors are listed here:
# https://www.rdkit.org/docs/GettingStartedInPython.html#list-of-available-descriptors

def get_parser():
    parser = argparse.ArgumentParser(
        usage="rdkit-descriptors.py < infile.tsv > outfile.tsv",
        description="Calculate physio-chemical properties from SMILES using RDKit.")
    parser.add_argument("smiles", nargs='?', default='smiles', help="Name of column with SMILES.")
    parser.add_argument("--cid", default='cid', help="Name of column with identifier for chemical. Default \"cid\".")
    parser.add_argument("--nproc", "-t", default=8, help="Number of parallel threads for computations.")
    parser.add_argument("--pdb", help="Write structures to PDB files in this directory.")
    parser.add_argument("--retries", default=10, type=int, help="Number of times to retry embedding conformers.")
    return parser

debug = False
if debug:
    # from src.utilities.io_utils import git_root
    # df = pd.read_table(git_root("results/5-chemicalFeatures/acceptors-props/pubchem.tsv"))
    args = get_parser().parse_args(["SMILES"])
    df = pd.read_table("~/Documents/phd/kittycat/chem/smiles.tsv")
else:
    args = get_parser().parse_args()
    df = pd.read_table(sys.stdin)

if args.pdb:
    os.makedirs(args.pdb, exist_ok=True)

# RDKit function names
funcs_core = ["NumValenceElectrons", "NumHeteroatoms",
              "MinPartialCharge", "MaxPartialCharge", "MinAbsPartialCharge", "MaxAbsPartialCharge",
              "MolWt", "NumHAcceptors", "NumHDonors", "NumRotatableBonds",
              # BertzCT is the complexity measure used by pubchem according to their glossary but seems a bit different.
              # There is a note about a fix done by RDKit so their value should be better.
              "BertzCT"]
funcs_MolDesc = ["CalcHallKierAlpha", "CalcLabuteASA",
                 "CalcNumRings", "CalcNumAmideBonds", "CalcNumBridgeheadAtoms", "CalcNumHeterocycles",
                 "CalcNumLipinskiHBA", "CalcNumLipinskiHBD", "CalcNumSpiroAtoms",
                 "CalcTPSA",
                 "CalcNumAliphaticRings", "CalcNumAliphaticCarbocycles", "CalcNumAliphaticHeterocycles",
                 "CalcNumAromaticRings", "CalcNumAromaticCarbocycles", "CalcNumAromaticHeterocycles",
                 "CalcNumSaturatedCarbocycles", "CalcNumSaturatedHeterocycles", "CalcNumSaturatedRings",
                 "CalcFractionCSP3"]
# requires a 3D conformer
funcs_3D = ["CalcEccentricity", "CalcAsphericity", "CalcInertialShapeFactor", "CalcNPR1", "CalcNPR2", "CalcSpherocityIndex", "CalcRadiusOfGyration"]
# descriptors that are not scalar. Currently not used.
funcs_array = ["BCUT2D"]

def compute_3D(m, row):
    if args.cid in row:
        desc = f"CID={row[args.cid]} "
    else:
        desc = ""
    desc += f"SMILES={row[args.smiles]}"
    AllChem.EmbedMolecule(m)
    for _ in range(args.retries):
        if m.GetNumConformers() > 0:
            return True
        AllChem.EmbedMolecule(m)
    sys.stderr.write(f"Failure: {desc}\n")
    return False

def pool_func(row):
    out = {}

    m = Chem.MolFromSmiles(row[1][args.smiles])
    m = Chem.AddHs(m)

    has3d = compute_3D(m, row[1])

    if has3d and args.pdb:
        id = row[1].get(args.cid, row[0])
        Chem.rdmolfiles.MolToPDBFile(m, os.path.join(args.pdb, str(id) + '.pdb'))

    for f in funcs_core:
        out[f] = getattr(Descriptors, f)(m)
    for f in funcs_MolDesc:
        colname = f
        if colname.startswith('Calc'): colname = colname[len('Calc'):]
        out[colname] = getattr(Descriptors.rdMolDescriptors, f)(m)
    for f in funcs_3D:
        colname = f
        if colname.startswith('Calc'): colname = colname[len('Calc'):]
        if has3d:
            out[colname] = getattr(Descriptors.rdMolDescriptors, f)(m)
        else:
            out[colname] = np.nan

    out["MolLogP"], out["MolMR"] = Descriptors.rdMolDescriptors.CalcCrippenDescriptors(m)
    out["HeavyAtomCount"] = m.GetNumHeavyAtoms()
    out["VanDerWaalsVolume"] = AllChem.ComputeMolVolume(m) if has3d else np.nan
    out["Charge"] = Chem.GetFormalCharge(m)
    out["SASA"] = rdFreeSASA.CalcSASA(m, rdFreeSASA.classifyAtoms(m)) if has3d else np.nan
    out["AtomStereoCount"] = rd.calcNumAtomStereo(m)
    out["BondStereoCount"] = rd.calcNumBondStereo(m)

    # from http://rdkit.org/docs/source/rdkit.Chem.Fragments.html
    for name, f in rd.numFragments.items():
        out[name] = f(m)

    for name in ['primary_alcohol', 'secondary_alcohol', 'tertiary_alcohol', 'coumarin', 'thiophenol', 'hydroxyamino']:
        out[name] = rd.get_num_substructs(m, rd.substructs[name])

    return out

if __name__ == '__main__':
    with Pool(args.nproc) as pool:
        rows = pool.map(pool_func, df.iterrows())

    out = pd.concat([df, pd.DataFrame(rows)], axis=1)
    # for backwards compatibility, keep names used when values were curated from pubchem
    out = out.rename(columns=dict(
        MolWt = "MolecularWeight",
        BertzCT = "Complexity",
        NumHDonors = "HBondDonorCount",
        NumHAcceptors = "HBondAcceptorCount",
        NumRotatableBonds = "RotatableBondCount"
    ))

    out.to_csv(sys.stdout, sep='\t', index=False)

