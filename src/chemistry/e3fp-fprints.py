#!/usr/bin/env python3
import sys, os
import argparse
from time import time
import pandas as pd
from e3fp.pipeline import fprints_from_smiles
from e3fp.fingerprint.db import FingerprintDatabase
from multiprocess import Pool, cpu_count
import logging
# avoid verbose logging
logging.basicConfig(level=logging.ERROR, force=True)

def pprint(*args):
    print(*args)
    sys.stdout.flush()

def get_parser():
    parser = argparse.ArgumentParser(description="Create E3FP fingerprints from smiles and add to new or existing database.")
    parser.add_argument("outfile", help="Filename of output, with (optional) extension .fpz.")
    parser.add_argument("-s", "--smiles", default="smiles", help="Name of column with SMILES. Default=\"smiles\".")
    parser.add_argument("-i", "--id", default="cid", help="Name of column with molecule identifiers. Default=\"cid\".")
    parser.add_argument("-t", "--threads", type=int, default=cpu_count(), help="Number of threads for multiprocess. Default=all.")
    parser.add_argument("-c", "--conformers", default=1, const=-1, action="store_const", help="Write all conformers for each id. Default=write only one.")
    parser.add_argument("--first", type=int, default=3, help="Max number of first conformers to generate fingerprints for.")
    parser.add_argument("--bits", type=int, default=2**16, help="Store fingerprints with this many bits. Must be a power of two.")
    return parser

debug = False
if debug:
    args = get_parser().parse_args([os.path.expanduser("~/Documents/phd/kittycat/chem/smiles-props/E3FP.fpz"), "-c"])
    df = pd.read_table("~/Documents/phd/kittycat/chem/smiles-props/pubchem.tsv")
else:
    args = get_parser().parse_args()
    df = pd.read_table(sys.stdin)

assert args.bits.bit_count() == 1, "--bits not a power of two"

df.set_index(args.id, inplace=True)

confgen_params = dict(first=args.first, num_conf=args.conformers)
fprint_params = dict(first=args.first, bits=args.bits)

if not args.outfile.endswith(".fpz"):
    args.outfile = args.outfile + ".fpz"
if os.path.isfile(args.outfile):
    db = FingerprintDatabase.load(args.outfile)
    completed = set([name.rsplit('_', maxsplit=1)[0] for name in db.fp_names])
    remaining = set(df.index) - completed
    pprint(f"# Read fingerprints. Completed: {len(completed)}. Remaining: {len(remaining)}.")
    df = df.loc[list(remaining)]
else:
    # level 5 is default. Adding new fingerprints will error if this isn't reflected here.
    db = FingerprintDatabase(level=5)

smiles = list(df[args.smiles])
ids = list(df.index)

def generate(s_n):
    s, n = s_n
    try:
        fps = fprints_from_smiles(s, n, confgen_params=confgen_params, fprint_params=fprint_params)
    except ValueError: # errors if the smiles is too simple, e.g. [H+]
        return None
    return fps

last = time()
successes = 0

with Pool(args.threads) as pool:
    for fps in pool.imap_unordered(generate, zip(smiles, ids)):
        if fps is not None:
            pprint(fps[0].name) # for progress tracking
            if args.conformers == 1:
                db.add_fingerprints(fps[0])
            else:
                db.add_fingerprints(fps)
            successes += 1
            current = time()
            # backup every minute and report progress
            if current - last > 60:
                pprint(f"# Completed: {successes}/{len(smiles)}")
                last = current
                db.savez(args.outfile)
                pprint(f"# Progress saved to {args.outfile}")

db.savez(args.outfile)
pprint(f"# Successes: {successes}/{len(smiles)}")

