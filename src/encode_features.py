#!/usr/bin/env python3
import argparse
import sys
import math
import numpy as np
import pandas as pd
import logging as log
from degnutil import input_output as io

def get_parser():
    parser = argparse.ArgumentParser(description="Encode features for ML, e.g. convert sequences to numeric encoding, remove redundant features, indicate train vs test set. Read train set from stdin, pipe /dev/null to use no train set.")

    parser.add_argument("test", nargs='*',
            help="""Filename(s) of test data with same format as infile or if multiple then they should have the same columns when combined. 
                    Multiple test files are combined in an all vs all fashion, so e.g. one file contains enzyme features and the other acceptor features. 
                    Train and test data files will be concatenated with a new column \"set\" with values \"train\" and \"test\".""")
    parser.add_argument("-i", "--ignore", dest="ignore", nargs="+", help="Name of columns that are not features so should be ignored and simply copied to output as identification.")
    parser.add_argument("--dna", nargs="+", help="Names of columns that should be encoded with alphabet ACGT. Not implemented.")
    parser.add_argument("--rna", nargs="+", help="Names of columns that should be encoded with alphabet ACGU. Not implemented.")
    parser.add_argument("--aa", nargs="+", help="Names of columns that should be encoded with a standard amino acid alphabet.")
    parser.add_argument("--aa-encoding", help="Path to encoding matrix for AAs, e.g. data/NCBI/blosum62Amb.tsv.")
    parser.add_argument("--cv", type=int, help="Do k-fold cross-validation, provide k. This will add a column \"set\" with integers in range [0;k-1].")
    parser.add_argument("-k", "--keep-redundant", action="store_true", help="Keep redundant features, i.e. features that never vary.")

    return parser


# functions

def cv_samples(nrow, k):
    """
    CV sampling without replacement.
    :param nrow: number of data points
    :param k: int. k-fold CV
    :return: list of integers in range [0;k-1] with length nrow
    """
    # create group indexes in range(k) for each datapoint (row) in df and shuffle the group indexes.
    group = np.asarray(list(range(k)) * math.ceil(nrow / k))[:nrow]
    np.random.shuffle(group)
    return group


def encode_seqs(seqs, encoding):
    """
    Encode string sequences with an encoding matrix, e.g. BLOSUM
    :param seqs: vector of strings
    :param encoding: pandas dataframe with alphabet as
    :return: np matrix with #row==len(seqs) and #col==len(seqs[0])*len(alphabet+1)
    """
    # make sure the gap character is the same.
    if "*" in encoding and "-" not in encoding and not any("*" in s for s in seqs):
        encoding = encoding.rename(columns={"*": "-"})
    if "-" in encoding and "*" not in encoding and not any("-" in s for s in seqs):
        encoding = encoding.rename(columns={"-": "*"})

    return np.asarray([np.asarray(encoding[list(s)]).flatten() for s in seqs])


def encode_seqs_df(seqs, encoding):
    """
    Same as encode_seqs but add column names to keep track of sequence code and sequence position.
    :param seqs: vector of strings. All strings has the same length otherwise the function will fail to make a matrix.
    :param encoding: matrix
    :return: pandas
    """
    # take len of encoding since the encoding matrix has to be square. If it is not, the extra columns should be ambiguity columns that e.g. represent multiple amino acids.
    alphabet = list(encoding)[:len(encoding)]
    seqlen = len(list(seqs)[0])
    fmt = "{:s}{:0" + str(len(str(seqlen))) + "d}"
    # we are using default flatten order="C" which is row-major, meaning we concatenate rows.
    # This means the alphabet along columns changes faster that the index for the position.
    names = [fmt.format(c, p) for p in range(seqlen) for c in alphabet]
    return pd.DataFrame(encode_seqs(seqs, encoding), columns=names)


def encode_seqs_df_col(df, col, encoding):
    """
    Same as encode_seqs_df except we find the seqs in a specific column col which is replaced by the encoded version.
    Naming is col + "_" + code + position
    :param df: pandas
    :param col: string name of column with sequences
    :param encoding: matrix
    :return: pandas
    """
    enc = encode_seqs_df(df[col], encoding).rename(columns=lambda n: col + "_" + n)
    # if filtering have happened, the df will have indexes that is not necessarily a simple range(0,?).
    # We could use pd.concat(..., ignore_index=True) but let's have the information pass through this function.
    enc.index = df.index
    return pd.concat([df, enc], axis=1, sort=False).drop(columns=col)


def identity_threshold_idx(train, test, threshold, split=None):
    """
    Get index of rows in train that satisfies a threshold of allowed identity with rows from test.
    :param train: pandas with a data point per row
    :param test: pandas with a data point per row
    :param threshold: float [0;1]
    :param split: optional list of names of columns which values should be split into each letter before doing identity comparisons
    :return: logical vector
    """
    train_mat = np.asarray(train.drop(columns=split))
    test_mat = np.asarray(test.drop(columns=split))

    if split is not None:
        for s in split:
            train_mat = np.hstack((train_mat, np.asarray([list(st) for st in train[s]])))
            test_mat = np.hstack((test_mat, np.asarray([list(st) for st in test[s]])))

    return np.asarray([np.mean(row == test_mat, axis=1).max() for row in train_mat]) <= threshold


def remove_redundant(df, ignore=None):
    """
    Find and remove columns where all values are identical across all given dfs
    :param df: pandas dataframe
    :param ignore: names of columns to ignore
    :return: list of strings names for columns that are redundant to use as features
    """
    df2 = df.drop(columns=ignore)
    redundant = df2.columns[np.all(df2 == df2.iloc[0, :], axis=0)]
    df.drop(columns=redundant, inplace=True)
    log.info("Removed {} out of {} ({:.2f}%) features that never varies".
             format(len(redundant), len(df2.columns), len(redundant) / len(df2.columns) * 100))
    return df


def encode_AAs(df, aa_encoding, aa_column_names):
    if aa_column_names is not None:
        for aa in aa_column_names:
            if aa in df:
                log.info(f"Encoding amino acids from {aa}.")
                df = encode_seqs_df_col(df, aa, aa_encoding)
    return df


def main(args):
    log.basicConfig(format="%(message)s", level=log.INFO)

    if args.aa is not None:
        if args.aa_encoding is None: raise NotImplementedError("Provide encoding matrix.")
        aa_encoding = pd.read_table(args.aa_encoding).select_dtypes(np.number)
    
    if args.infile is None:
        df = None
    else:
        try: df = encode_AAs(pd.read_table(args.infile), aa_encoding, args.aa)
        except pd.errors.EmptyDataError: df = None

    if len(args.test) > 0:
        if args.cv is not None: raise NotImplementedError("Test file not used with CV.")
        test = encode_AAs(pd.read_table(args.test[0]), aa_encoding, args.aa)
        if len(args.test) > 1:
            sizes = [len(test)]
            for f in args.test[1:]:
                test_right = encode_AAs(pd.read_table(f), aa_encoding, args.aa)
                sizes.append(len(test_right))
                test = test.merge(test_right, how="cross")
            assert len(test) == np.prod(sizes)
            log.info(f"{len(sizes)} files cross combined resulting in {'*'.join(map(str, sizes))} = {len(test)} data points.")

        if df is not None:
            df['set'] = "train"
            test['set'] = "test"
            df = pd.concat([df, test])
            args.ignore += ['set']
        else:
            df = test

    if args.cv is not None:
        df['set'] = cv_samples(len(df), args.cv)
        args.ignore += ['set']

    # for now drop string features, but in the future they should maybe be implemented as 1-hot automatically
    feature_names = np.setdiff1d(df.columns, args.ignore)
    feature_dtypes = df[feature_names].dtypes
    nonum = feature_names[(feature_dtypes != int) & (feature_dtypes != float)]
    if len(nonum) > 0:
        log.info("Dropping non-number feature(s): " + ','.join(nonum))
        df.drop(columns=nonum, inplace=True)

    if not args.keep_redundant:
        remove_redundant(df, args.ignore)

    df.to_csv(sys.stdout, sep='\t', index=False)



if __name__ == '__main__':
    args = get_parser().parse_args()
    args.infile = sys.stdin if io.is_reading_from_pipe() else None
    main(args)

debug = False
if debug:
    args = get_parser().parse_args(
        "~/biosustain/gt/features/rdkit-desc_muscle-ntern-hmm_median_notgtpred.tsv -i reaction enzyme acceptor cid source rate --aa seq --aa-encoding ~/biosustain/data/ncbi/match.tsv".split())
    args.infile = "~/biosustain/gt/features/rdkit-desc_muscle-ntern-hmm_median_gtpred.tsv"
    main(args)

