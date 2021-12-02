#!/usr/bin/env python3
import argparse
import sys
import numpy as np
import pandas as pd
from rdkit import Chem
import src.chemistry.chem_utils.rdkit_utils as u_rdkit
import re
import logging


def get_parser():
    parser = argparse.ArgumentParser(
            usage="mol2pdf.py identifiers.tsv [OPTIONS]", 
            description="Read a table of identifiers, draw structures.")
    parser.add_argument('infile', nargs='?', default=sys.stdin, help='Name of infile. Default=read from stdin, in which case the -o/--out is required.')
    parser.add_argument('-o', '--out', help='Name of outfile. Default=infile with extension replaced by ".pdf".')
    parser.add_argument('-S', '--smiles', help='Name of column with SMILES to use. Default is case-insensitive match to "SMILES".')
    parser.add_argument('-C', '--cas', help='Name of column with CAS to use. Default is case-insensitive match to "CAS".')
    parser.add_argument('-I', '--inchi', help='Name of column with INCHI to use. Default is case-insensitive match to "INCHI".') 
    parser.add_argument('-K', '--inchikey', help='Name of column with INCHI key to use. Default is case-insensitive letter-only match to "INCHIKEY".')
    parser.add_argument('--names', nargs="+", help='Name of columns to use as names in the plot. Default is case-insensitive match of "NAME", "ID" and using CAS if found. Otherwise fallback on other identifiers.')
    return parser

def get_args(args_list=None):
    args = get_parser().parse_args(args_list)
    if args.out is None:
        if args.infile == sys.stdin:
            raise argparse.ArgumentTypeError("Outfile is required when using stdin.")
        args.out = args.infile.rsplit('.', 1)[0] + '.pdf'
    return args


debug = True
if debug:
    args = get_args(["~/biosustain/gt/misc/L6100_polyphenolicLib.tsv"])
else:
    args = get_args()


df = pd.read_table(args.infile)
cols_upper = [re.sub('[^A-Z]', '', c.upper()) for c in df.columns]

identifiers = {}
for key in ["id", "name", "smiles", "cas", "inchi", "inchikey"]:
    explicit = getattr(args, key, None)
    if explicit is not None:
        identifiers[key] = df[explicit]
    else:
        try: identifiers[key] = df.iloc[:, cols_upper.index(key.upper())]
        except ValueError: pass

# find names and identifiers to show
N = len(df)
names = [{} for _ in range(N)]
for key in ["id", "name", "cas"]:
    if key in identifiers:
        for i in range(N):
            names[i][key] = identifiers[key][i]

# fallback to other identifiers
for i in range(N):
    if len(names[i]) == 0:
        for key in ["inchikey", "smiles", "inchi"]:
            try: names[i][key] = identifiers[key][i]
            except KeyError: pass
            else: break


# fallback to a simple count
for i in range(N):
    if len(names[i]) == 0:
        names[i]["index"] = i

names = ['\n'.join(f'{k}={v}' for k, v in d.items()) for d in names]

mols = None

if "smiles" in identifiers:
    mols = u_rdkit.smiles2molecules(identifiers["smiles"], names=names)
    u_rdkit.molecules2pdf(args.out, mols)

if "cas" in identifiers:
    mols_cas = u_rdkit.identifier2molecules(identifiers["cas"], names=names)
    mols_cas = [m for m in mols_cas if m is not None]
    u_rdkit.molecules2pdf(args.out.rsplit('.', 1)[0] + '_cas.pdf', mols_cas)

if "inchi" in identifiers:
    logging.warning("INCHI not used to make molecules. Not implemented yet.")

if "inchikey" in identifiers:
    logging.warning("INCHI key not used to make molecules. Not implemented yet.")



