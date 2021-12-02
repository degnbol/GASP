#!/usr/bin/env python3
import argparse
import sys
import numpy as np
import pandas as pd
import logging as log
from sklearn import svm
from sklearn.metrics import roc_auc_score as AUC


def get_parser():
    parser = argparse.ArgumentParser(
        description="Train a support vector machine and perform prediction.")

    parser.add_argument("-i", "--ignore", default=[], dest="ignore", nargs="+", help="Name of columns that are neither feature nor class so are simply copied to output for identification.")
    parser.add_argument("-c", "--class", dest="Class", help="Name of column with class labels in training data. Default=class.", default="class")
    parser.add_argument("-s", "--set", help="Name of column with indication of set membership, e.g. \"train\" and \"test\". Default=set.", default="set")

    return parser


# functions


def main(args):
    log.basicConfig(format="%(message)s", level=log.INFO)

    # reading

    df = pd.read_table(args.infile)
    trainset = df[df[args.set] == "train"]
    testset = df[df[args.set] == "test"]

    clf = svm.SVC(probability=True)
    clf.fit(trainset.drop(columns=args.ignore + [args.Class, args.set]), trainset[args.Class])
    pred = clf.predict_proba(testset.drop(columns=args.ignore + [args.Class, args.set]))

    print(AUC(testset[args.Class], pred[:,1]))



if __name__ == '__main__':
    args = get_parser().parse_args()
    args.infile = sys.stdin
    main(args)


debug = False
if debug:
    args = get_parser().parse_args("-i enzyme acceptor source cid rate -c reaction".split())
    args.infile = "~/biosustain/gt/randomforest/rdkit-desc_muscle-ntern-hmm_median_gtpred/traintest.tsv"
    main(args)
