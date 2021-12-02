#!/usr/bin/env python3
import sys
import argparse
import pandas as pd
from e3fp.pipeline import fprints_from_smiles
from e3fp.fingerprint.db import FingerprintDatabase
from python_utilities.parallel import Parallelizer
import logging
# avoid verbose logging
logging.basicConfig(level=logging.WARNING, force=True)

def get_parser():
    parser = argparse.ArgumentParser(description="Create E3FP fingerprints from smiles.")
    parser.add_argument("outfile", help="Filename of output, e.g. ending in .fpz. We can't use stdout since e3fp code uses the file ending in the savez code.")
    parser.add_argument("-s", "--smiles", default="smiles", help="Name of column with SMILES. Default=\"smiles\".")
    parser.add_argument("-i", "--id", default="cid", help="Name of column with molecule identifiers. Default=\"cid\".")
    return parser

debug = False
if debug:
    df = pd.read_table("/gt/acceptors/acceptor_smiles.tsv")
    smiles = df.smiles
    ids = df.cid
    outfile = None
else:
    args = get_parser().parse_args()
    df = pd.read_table(sys.stdin)
    smiles = df[args.smiles]
    ids = df[args.id]
    outfile = args.outfile

parallel = Parallelizer("processes")
fprints = parallel.run(fprints_from_smiles, ((s, n) for s, n in zip(smiles, ids)))
# level is number of iterations performed before stopping generation of fingerprints. Just have to match.
# Default is 5. -1 means until termination.
db = FingerprintDatabase(name="acceptors", level=5)
db.add_fingerprints([f for fs in fprints for f in fs[0]])

db.savez(outfile)
