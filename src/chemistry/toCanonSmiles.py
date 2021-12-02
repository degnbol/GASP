#!/usr/bin/env python3
import sys
import pandas as pd
from rdkit import Chem
import argparse

def get_args():
    parser = argparse.ArgumentParser(usage="src/chemistry/toCanonSmiles.py [options] < infile > outfile",
        description="Convert e.g. smiles from pubchem to canonical smiles made with RDKit. In future could also do inchi etc.")
    parser.add_argument("-H", "--header", action="store_true", help="Does infile have a header? Default is an id on each line.")
    parser.add_argument("-i", "--col", default="smiles", help="Name of column with input used if -H/--header.")
    parser.add_argument("-d", "--sep", default='\t', help="Delimiter if input is a table.")
    parser.add_argument("-o", "--newcol", help="Print infile with the result in a new column. Default is only printing the result.")
    return parser.parse_args()

args = get_args()
if args.header:
    df = pd.read_table(sys.stdin, sep=args.sep)
    ids = df[args.col]
else:
    with sys.stdin as infile:
        ids = [l.strip() for l in infile]

canonSmiles = [Chem.CanonSmiles(i) for i in ids]

if args.newcol is None:
    print('\n'.join(canonSmiles))
else:
    if not args.header: raise NotImplementedError
    df[args.newcol] = canonSmiles
    df.to_csv(sys.stdout, sep=args.sep, index=False)
