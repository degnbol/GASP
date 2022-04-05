#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
from src.chemistry.chem_utils.rdkit_utils import *
from src.utilities.io_utils import git_root

# Generates negatives from positives by applying small modifications according to the rules listed below.
# A negative is a chemical that is reactive with non of the GT1 enzymes. 
# The positives are chemicals with experimental data validating that they are reactive with at least one GT1 enzyme.

# Rules for generating negative acceptors:
# 1) if the molecule presents a -NH group (R-NH2, or R-NH-R'), ignore it.
# 2) In other molecules, replace each OH-and -SH groups by -H (best negative analogue), and
# 3) ideally a random mix of -H, -F, OMe (OCH3), -Cl, OAc (OCOCH3).
# As an example, ethanol should become ethane, 1-fluoro ethane, methoxyethane, Ethyl acetate.

# USAGE: Run the script ./gen_negatives.py or interactively. SMILES from positives must be present at the assumed path.

os.chdir(git_root("results/3-generateNegatives"))
try: os.mkdir("diagrams") # for output .pdfs, .svgs, .pngs
except OSError: pass

# Validation; does simple compounds look as expected?
ethanol = Chem.MolFromSmiles('CCO')
ethane, = Chem.ReplaceSubstructs(ethanol, hydroxy, hydrogen)
fluoro_ethane, = Chem.ReplaceSubstructs(ethanol, hydroxy, flour)
methoxyethane, = Chem.ReplaceSubstructs(ethanol, hydroxy, OMe)
ethyl_acetate, = Chem.ReplaceSubstructs(ethanol, hydroxy, OAc)
# render_pdf("diagrams/ethanol.pdf", [ethanol, ethane, fluoro_ethane, methoxyethane, ethyl_acetate], legends=["ethanol", "ethane", "1-fluoro ethane", "methoxyethane", "ethyl acetate"])

# SMILES from chemicals that are experimentally observed reactive for at least one GT1
positives = pd.read_table("positives.tsv")

mols = smiles2molecules(*(positives[k] for k in ["smiles", "cid"]))
mols = mols[~has_frag(mols, "NH0", "NH1", "NH2")] # filter out NH
# Not sure if all nitrogen containing molecules should be removed but that decision is not necessary to make as there are non left:
assert not has_atom(mols, 'N').any(), "things changed"

mols = COOH2COOMe(mols)
mols_H = replace_sub(replace_sub(mols, "C[OH]", "C"), "c[OH]", "c")

assert not has_frag(mols_H, "Al_OH", "Ar_OH").any(), "Not all R-OH has been removed."
# render_pdf("diagrams/negatives_H.pdf", mols_H)

# random replacement versions
assert not has_hydroxy(mols_H).any(), "The current version of randomly_replace_hydroxy will replace ANY OH, not just Al_OH and Ar_Oh"

mols = mols[np.argsort([int(c) for c in get_names(mols)])]
mols_rnd = randomly_replace_OH_SH(mols)
mols_rnd = [Molecule(m) for m in mols_rnd]
suffixes = ["H", "F", "Cl", "OMe", "OAc"]
for i, mol in enumerate(mols_rnd): mol.name += "_" + suffixes[i % len(suffixes)]
molecules2pdf("negatives_rnd.pdf", mols_rnd)

with open("negatives.smiles", 'w') as fh:
    fh.write('\n'.join({m.to_smiles() for m in mols_rnd}))

