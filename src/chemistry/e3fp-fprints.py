#!/usr/bin/env python3
import sys
import argparse
import pandas as pd
from e3fp.pipeline import fprints_from_smiles
from e3fp.fingerprint.db import FingerprintDatabase
# from python_utilities.parallel import Parallelizer
import multiprocess as mp
import logging
# avoid verbose logging
logging.basicConfig(level=logging.WARNING, force=True)

def get_parser():
    parser = argparse.ArgumentParser(description="Create E3FP fingerprints from smiles.")
    parser.add_argument("outfile", help="Filename of output, e.g. ending in .fpz. We can't use stdout since e3fp code uses the file ending in the savez code.")
    parser.add_argument("-s", "--smiles", default="smiles", help="Name of column with SMILES. Default=\"smiles\".")
    parser.add_argument("-i", "--id", default="cid", help="Name of column with molecule identifiers. Default=\"cid\".")
    parser.add_argument("-t", "--threads", type=int, default=mp.cpu_count(), help="Number of threads for multiprocess. Default=all.")
    return parser

debug = False
if debug:
    from src.utilities.io_utils import git_root
    df = pd.read_table(git_root("results/5-chemicalFeatures/20220215-props/pubchem.tsv"))
    args = get_parser().parse_args(["E3FP.tmp.fpz"])
else:
    args = get_parser().parse_args()
    df = pd.read_table(sys.stdin)

smiles = df[args.smiles]
ids = df[args.id]

def progress(s, n):
    fps = fprints_from_smiles(s, n)
    print(n)
    return fps

with mp.Pool(args.threads) as pool:
    fprints = pool.starmap(progress, zip(smiles, ids))

# level is number of iterations performed before stopping generation of fingerprints. Just have to match.
# Default is 5. -1 means until termination.
db = FingerprintDatabase(name="acceptors", level=5)
# NOTE: only adding one conformer
db.add_fingerprints([fs[0] for fs in fprints])
db.savez(args.outfile)

