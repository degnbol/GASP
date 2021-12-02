#!/usr/bin/env python3
import numpy as np
from os.path import expanduser
from e3fp.fingerprint.db import FingerprintDatabase as DB
from e3fp.fingerprint.metrics.fprint_metrics import distance
from sklearn.manifold import MDS
import argparse

def get_parser():
    parser = argparse.ArgumentParser(usage="%(prog)s infile.fpz [options] > outfile.tsv",
             description="Convert numpy compressed (npz format) fingerprints made with E3FP to points in an k-dimensional space. "
                         "Fingerprints are converted by calculating all pairwise Euclidean distances followed by metric MultiDimensional Scaling.")
    parser.add_argument("infile", help="Filename for chemical fingerprints made with E3FP, e.g. ending in fpz.")
    parser.add_argument("-k", "--dimensions", default=2, type=int, help="Number of dimensions to reduce the fingerprint bitvectors to.")
    parser.add_argument("-c", "--conformers", action="store_true", help="Set distance to zero between multiple conformers made for the same id.")
    return parser

debug = False
if debug:
    args = get_parser().parse_args("~/biosustain/gt/acceptors/acceptors.fpz -ck 12".split(' '))
    main(args)

def main(args):
    db = DB.load(expanduser(args.infile))
    names = np.asarray([n.rsplit('_', 1)[0] for n in db.fp_names])
    distances = np.asarray([[distance(fp1, fp2) for fp2 in db] for fp1 in db])

    if args.conformers:
        # modify distances so chemicals with multiple fingerprints (multiple conformers) have zero distance between the different fingerprints
        for i, n in enumerate(names):
            distances[i, names == n] = 0

    mds = MDS(n_components=args.dimensions, dissimilarity="precomputed")
    mds.fit(distances)
    embeddings = {n: np.mean(mds.embedding_[names == n, :], axis=0) for n in set(names)}

    # write result
    decimals = len(str(args.dimensions))
    feature_fmt = "MDS_{:0%dd}" % decimals
    print("id\t" + '\t'.join([feature_fmt.format(i+1) for i in range(args.dimensions)]))
    for n, vals in embeddings.items():
        print(n, *vals, sep='\t')


if __name__ == '__main__':
    main(get_parser().parse_args())
