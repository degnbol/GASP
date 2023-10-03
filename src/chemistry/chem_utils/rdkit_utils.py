#!/usr/bin/env python3
import os, sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import Fragments, Draw
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF
import cirpy


class Molecule(Chem.Mol):
    """
    Extension of RDKit molecule that checks identity and uniqueness on SMILES alone.
    """

    def __init__(self, smiles_or_molecule, name=None):
        if type(smiles_or_molecule) == str: smiles_or_molecule = Chem.MolFromSmiles(smiles_or_molecule)
        super().__init__(smiles_or_molecule)
        if name is not None: self.name = name

    def to_smiles(self):
        return Chem.MolToSmiles(self)

    @property
    def name(self):
        return self.GetProp("_Name")

    @name.setter
    def name(self, value):
        self.SetProp("_Name", str(value))

    def __eq__(self, other):
        if type(other) != Molecule: return False
        return self.to_smiles() == other.to_smiles()

    def __hash__(self):
        return hash("Mol:"+self.to_smiles())





### CONSTANTS

hydroxy = Chem.MolFromSmarts("[OH]")
thiol = Chem.MolFromSmarts("[SH]")
Al_OH = Chem.MolFromSmarts("C[OH]")
Ar_OH = Chem.MolFromSmarts("c[OH]")
carboxyl = Chem.MolFromSmarts("C(=O)[OH]")
carboxylate = Chem.MolFromSmarts("C(=O)[O-]")
hydrogen = Chem.MolFromSmarts("[H]")
flour = Chem.MolFromSmarts("F")
chlorine = Chem.MolFromSmarts("Cl")
OMe = Chem.MolFromSmarts("O[CH3]")
OAc = Chem.MolFromSmarts("OC(=O)[CH3]")
COOMe = Chem.MolFromSmarts("C(=O)O[CH3]")

substructs = dict(OMe=OMe, OAc=OAc, COOMe=COOMe)

MOLDIR = os.path.join(os.path.dirname(__file__), '..', 'mols')
for fname in os.listdir(MOLDIR):
    if fname.endswith(".mol"):
        name = fname[:-len(".mol")]
        substructs[name] = Chem.MolFromMolFile(f"{MOLDIR}/{fname}")


numFragments = {k[len("fr_"):]: f for k, f in Fragments.fns}


def smiles2molecules(smiles, names=None):
    molecules = [Chem.MolFromSmiles(s) for s in smiles]
    if names is not None: set_names(molecules, names)
    return np.asarray(molecules)


def canonical(molecules):
    """
    Convert smiles within molecules to canonical SMILES.
    """
    out = [Molecule(Chem.MolFromSmiles(Chem.MolToSmiles(m))) for m in molecules]
    try: names = get_names(molecules)
    except KeyError: pass
    else: set_names(out, names)
    return np.asarray(out)
    

def identifier2molecules(queries, names=None):
    """
    Convert CAS or other identifiers to molecules using CIRpy (Chemical Identifier Resolver).
    :param queries: input queries, e.g. CAS numbers.
    :param names: list of names to give the molecules if found.
    :return: list of molecules and None when there was no match.
    """
    molecules = []
    for i, q in enumerate(queries):
        print(f"{i}/{len(queries)}")
        s = cirpy.resolve(q, "smiles")
        if s is None: molecules.append(None)
        else:
            molecules.append(Chem.MolFromSmiles(s))

    if names is not None: set_names(molecules, names)
    return molecules


def calcNumAtomStereo(mol):
    return sum(s.type.name.startswith("Atom") for s in Chem.FindPotentialStereo(mol))
def calcNumBondStereo(mol):
    return sum(s.type.name.startswith("Bond") for s in Chem.FindPotentialStereo(mol))


def _has_frag(molecules, fragment):
    return np.asarray([numFragments[fragment](m) > 0 for m in molecules])
def has_frag(molecules, *fragments):
    out = np.zeros(len(molecules), dtype=bool)
    for fragment in fragments: out |= _has_frag(molecules, fragment)
    return out

def get_unique_atoms(molecule):
    return {a.GetSymbol() for a in molecule.GetAtoms()}

def has_atom(molecules, atom):
    return np.asarray([atom in get_unique_atoms(m) for m in molecules])

def get_num_substructs(molecule, substructure):
    return len(molecule.GetSubstructMatches(substructure))

def has_sub(molecules, query):
    if type(query) == str: query = Chem.MolFromSmarts(query)
    return np.asarray([molecule.HasSubstructMatch(query) for molecule in molecules])

def has_any_sub(molecules, *queries):
    out = np.zeros(len(molecules), dtype=bool)
    for sub in queries: out |= has_sub(molecules, sub)
    return out

def has_hydroxy(molecules):
    return has_sub(molecules, hydroxy)

def has_carboxylate(molecules):
    return has_sub(molecules, carboxylate)

def get_name(molecule):
    return molecule.GetProp("_Name")
def get_names(molecules):
    return [get_name(molecule) for molecule in molecules]

def get_smiles(molecules):
    return [Chem.MolToSmiles(molecule) for molecule in molecules]

def set_name(molecule, name):
    molecule.SetProp("_Name", str(name))
def set_names(molecules, names):
    for molecule, name in zip(molecules, names):
        if molecule is not None:
            set_name(molecule, name)

def _replace_sub(molecule, query, replacement):
    return Chem.ReplaceSubstructs(molecule, query, replacement, replaceAll=True)[0]
def replace_sub(molecules, query, replacement):
    if type(query) == str: query = Chem.MolFromSmarts(query)
    if type(replacement) == str: replacement = Chem.MolFromSmarts(replacement)
    return np.asarray([_replace_sub(molecule, query, replacement) for molecule in molecules])

def COO2COOH(molecules):
    return replace_sub(molecules, carboxylate, carboxyl)

def COOH2COOMe(molecules):
    return replace_sub(molecules, carboxyl, COOMe)

def remove_Hs(molecules):
    return np.asarray([Chem.RemoveAllHs(molecule) for molecule in molecules])

def remove_sub(molecules, query):
    if type(query) == str: query = Chem.MolFromSmarts(query)
    return np.asarray([Chem.DeleteSubstructs(molecule, query) for molecule in molecules])

def remove_any_sub(molecules, *queries):
    for query in queries: molecules = remove_sub(molecules, query)
    return molecules

def _remove_any_sub(molecule, *queries):
    for query in queries: molecule = Chem.DeleteSubstructs(molecule, query)
    return molecule

def replace_random(molecule, query, replacement):
    return Chem.RemoveAllHs(np.random.choice(Chem.ReplaceSubstructs(molecule, query, replacement)))

def replace_any_random(molecule, queries, replacement):
    """
    Replace one of the instances among queries in a molecule.
    If non of the queries are found in the molecule it will be returned unchanged.
    """
    queries = [q for q in queries if molecule.HasSubstructMatch(q)]
    if len(queries) == 0: return molecule
    choices = [c for q in queries for c in Chem.ReplaceSubstructs(molecule, q, replacement)]
    return Chem.RemoveAllHs(np.random.choice(choices))

def randomly_replace(molecule, query, replacements):
    """
    If query is not found anywhere in the molecule, it will be returned unchanged.
    :param molecule:
    :param query:
    :param replacements:
    :return:
    """
    while molecule.HasSubstructMatch(query):
        # pick a random replacement
        # replace in one of the possible locations
        molecule = replace_random(molecule, query, np.random.choice(replacements))
    return molecule

def randomly_replace_remove(molecule, query, replacement):
    """
    Make versions of the given molecule where all queries are removed,
    except for one randomly selected query that will be replaced with "replacement".
    :param molecule: molecule to make modified copies of
    :param query: molecule substructure to replace one of and remove the rest
    :param replacement: molecule
    :return: new molecule
    """
    return Chem.DeleteSubstructs(replace_random(molecule, query, replacement), query)

def randomly_replace_remove_any(molecule, queries, replacement):
    """
    Make versions of the given molecule where all queries are removed,
    except for one randomly selected query that will be replaced with "replacement".
    :param molecule: molecule to make modified copies of
    :param queries: molecule substructures to replace one of and remove the rest
    :param replacement: molecule
    :return: new molecule
    """
    return _remove_any_sub(replace_any_random(molecule, queries, replacement), *queries)


def randomly_replace_OH_SH(molecules):
    """
    Take any molecule with OH, then remove any OH and any SH except a randomly 
    selected one that will be replaced by one of the replacement choices.
    WARNING: All OH, so it will also be OH that are not Al_OH or Ar_OH, if present.
    :param molecules: vector of molecules
    :return: randomly substituted versions of the given molecules (potentially with duplicates)
    """
    queries = [hydroxy, thiol]
    molecules = np.asarray(molecules)[has_hydroxy(molecules)]
    replacements = [hydrogen, flour, chlorine, OMe, OAc]
    return np.asarray([randomly_replace_remove_any(m, queries, r) for m in molecules for r in replacements])


def molecules2pdf(path, molecules, legends=None, resolution=400):
    path = os.path.expanduser(path)
    molecules = [Chem.RemoveAllHs(m) for m in molecules]
    if legends is None: legends = get_names(molecules)
    legends = [str(l) for l in legends]
    svg_string = Draw.MolsToGridImage(molecules, molsPerRow=min(5, len(molecules)), subImgSize=(resolution, resolution), legends=legends, useSVG=True)
    with open(path, 'w') as outfile: outfile.write(svg_string)
    renderPDF.drawToFile(svg2rlg(path), path)



