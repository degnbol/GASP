#!/usr/bin/env python3
import argparse
import sys
import math
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score as AUC
import joblib
import logging as log


def get_parser():
    parser = argparse.ArgumentParser(
        description="Train a Random Forest classifier with a train and test set and make prediction on the test set. "
                    "Optionally use CV to split infile into train and test. "
                    "Reads a tab-separated table header from stdin and writes prediction table to stdout. "
                    "Columns: id, features, set, and class. Class is 1 for a positive and 0 for a negative. set indicates train or test. Create infile with encode_features.py.")

    parser.add_argument("-i", "--id", dest="ids", nargs="+", help="Name of columns that are neither feature nor class so should be identifier columns. They are simply copied to output for identification.")
    parser.add_argument("-F", "--features", nargs="+", help="Name of columns that are features to be used.")
    parser.add_argument("-t", "--identity", dest="threshold", type=float, default=1., help="Threshold for train vs test identity.")
    parser.add_argument("-o", "--out",
                        help="The model can be saved to this file (e.g. randomforest.joblib). By default the model is not saved.")
    parser.add_argument("-m", "--model",
                        help="Read model (e.g. randomforest.joblib) from this arg and test instead of train from stdin.")
    parser.add_argument("-n", "--estimators", dest="n_estimators", type=int, default=100, help="Number of trees in the forest.")
    parser.add_argument("--crit", dest="criterion", default="gini", help="Randomforest criterion for split quality.", choices=["gini", "entropy"])
    # uppercase Class since class is a reserved word
    parser.add_argument("-c", "--class", dest="Class", help="Name of column with class labels in training data. Default=class.", default="class")
    parser.add_argument("-e", "--eval", default="class",
                        help="Name of column with evaluation value in test data. Can be class like infile but floats are also supported. Default=class.")
    parser.add_argument("--eval-thres", dest="eval_threshold", type=float, help="Set this threshold to also turn a float eval from test set into two classes.")
    parser.add_argument("--log", help="Logging is sent to stderr. Set this flag to a filename to also write it to a file.")
    parser.add_argument("--importance", help="Write feature importance and stddev to this file.")
    parser.add_argument("--cv", type=int, help="Do k-fold cross-validation, provide k.")
    parser.add_argument("-s", "--set", help="Name of column with set indication. Can have comma separated values. The test set is each unique name found. Default=\"set\".", default="set")

    return parser


# functions

def identity_threshold_idx(train, test, threshold):
    """
    Get index of rows in train that satisfies a threshold of allowed identity with rows from test.
    :param train: pandas with a data point per row
    :param test: pandas with a data point per row
    :param threshold: float [0;1]
    :return: logical vector
    """
    train_mat = np.asarray(train)
    test_mat = np.asarray(test)
    return np.asarray([np.mean(row == test_mat, axis=1).max() for row in train_mat]) <= threshold


def remove_redundant(*dfs, ignore=None):
    """
    Find and remove columns where all values are identical across all given dfs
    :param dfs: pandas dataframes
    :param ignore: names of columns to ignore
    :return: list of strings names for columns that are redundant to use as features
    """
    # skip None
    dfs = [df for df in dfs if df is not None]
    cat = pd.concat(dfs).drop(columns=ignore)
    redundant = cat.columns[np.all(cat == cat.iloc[0, :], axis=0)]
    if len(redundant) > 0:
        log.info("Removed {} out of {} ({:.2f}%) features that never varies".
                 format(len(redundant), len(cat.columns), len(redundant) / len(cat.columns) * 100))
        return [df.drop(columns=redundant) for df in dfs]
    return dfs


def feature_prep(train, test, threshold=1.0, ignore=None):
    """
    Prepare features for random forest. This means we make sure to have less similarity between train and test than a threshold,
    we encode AA sequences with a given aa_encoding matrix, we drop redundant features that never vary.
    :param train: pd dataframe
    :param test: pd dataframe
    :param threshold: threshold that identity must be lower than
    :param ignore: columns to ignore that are not features
    :return: modified train and test set ready for random forest
    """
    if threshold < 1:
        trainidx = identity_threshold_idx(train.drop(columns=ignore, errors='ignore'),
                                          test.drop(columns=ignore, errors='ignore'), threshold)

        train = train.loc[trainidx, :]
        log.info("{} out of {} ({:.2f}%) training data points discarded due to identity > {}% with a test data point".
                 format(sum(~trainidx), len(trainidx), (1-trainidx.mean())*100, threshold*100))
    assert len(train) > 0, "All training data removed"

    # for now drop string features, but in the future they should maybe be implemented as 1-hot automatically
    train_feature_names = np.setdiff1d(train.columns, ignore)
    test_feature_names = np.setdiff1d(test.columns, ignore)
    train_feature_dtypes = train[train_feature_names].dtypes
    test_feature_dtypes = test[test_feature_names].dtypes
    train_nonum = train_feature_names[(train_feature_dtypes != int) & (train_feature_dtypes != float)]
    test_nonum = test_feature_names[(test_feature_dtypes != int) & (test_feature_dtypes != float)]
    assert np.all(train_nonum == test_nonum)
    if len(train_nonum) > 0:
        log.info("Dropping non-number feature(s): " + ','.join(train_nonum))
        train.drop(columns=train_nonum, inplace=True)
        test.drop(columns=test_nonum, inplace=True)

    train, test = remove_redundant(train, test, ignore=ignore)
    return train, test


def training(train, Class, n_estimators, criterion):
    """
    Perform training of random forest classifier with proper features given.
    :param train: training features.
    :param Class: name of column with class label.
    :param n_estimators:
    :param criterion:
    :return: trained RandomForestClassifier
    """
    clf = RandomForestClassifier(n_estimators=n_estimators, n_jobs=-1, criterion=criterion)
    clf.fit(train.drop(columns=Class), train[Class])
    return clf


def testing(clf, test, ignore):
    predictions = clf.predict_proba(test.drop(columns=ignore, errors='ignore'))[:, clf.classes_ == 1].squeeze()
    # explicitly call copy since you get a view otherwise and the assignment after fails
    test["pred"] = predictions
    return test


def traintest(train, test, Class, threshold, n_estimators, criterion, ignore=None):
    if ignore is None: ignore = []
    elif type(ignore) is not list: ignore = [ignore]
    train = train.drop(columns=ignore, errors='ignore')
    train, test = feature_prep(train, test, threshold, [Class] + ignore)
    log.info("Train random forest model.")
    clf = training(train, Class, n_estimators, criterion)
    if len(test) == 0: return None, clf
    log.info("Make prediction on the test set.")
    test = testing(clf, test, ignore + [Class])
    return test, clf


def cv_samples(df, k):
    """
    Yield random train and test samples of rows from df k times.
    :param df: pandas dataframe of features
    :param k: int. k-fold CV
    :return: train and test set for a random sample without replacement from df. len(df)/k rows are yielded
    """
    # create group indexes in range(k) for each datapoint (row) in df and shuffle the group indexes.
    group = np.asarray(list(range(k)) * math.ceil(len(df) / k))[:len(df)]
    np.random.shuffle(group)
    for g in range(k):
        yield df.iloc[group != g, :], df.iloc[group == g, :]


def PCC(x, y):
    return np.corrcoef(x, y)[0,1]


def print_performance(df, class_name, eval_name, pred_name, eval_threshold=None):
    Class = df.get(class_name)
    Eval = df.get(eval_name)
    Pred = df.get(pred_name)

    if Class is not None:
        log.info("PCC={:.3f} to {}".format(PCC(Pred, Class), class_name))
        class_auc = AUC(Class, Pred)
        log.info("AUC={:.3f} for {}/{} {} values".format(class_auc, sum(Class), len(Class), class_name))


    if Eval is not None:
        Pred, Eval = Pred[~np.isnan(Eval)], Eval[~np.isnan(Eval)]

        log.info("PCC={:.3f} to {}".format(PCC(Pred, Eval), eval_name))

        if Eval.dtype == float:
            response = Pred > 0.5  # majority vote
            resp_auc = AUC(response, Eval)
            log.info("AUC={:.3f} for {}/{} positive responses".format(resp_auc, sum(response), len(response)))
            if eval_threshold is not None:
                thres_eval = Eval > eval_threshold
                thres_auc = AUC(thres_eval, Pred)
                log.info("AUC={:.3f} for {}/{} {} values > {}".format(thres_auc, sum(thres_eval), len(thres_eval), eval_name, eval_threshold))
        else:
            log.info("AUC={:.3f}".format(AUC(Eval, Pred)))




def main(args):
    if args.log is None: handlers = None
    else: handlers = [log.StreamHandler(), log.FileHandler(args.log)]
    log.basicConfig(format="%(message)s", level=log.INFO, handlers=handlers)

    if args.set != "set":
        log.info("Using column {} as set".format(args.set))

    # reading

    infile = pd.read_table(args.infile)
    if infile.isna().any().any():
        log.warning("NA found in infile in column(s): " + ",".join(infile.columns[infile.isna().any()]))
    if args.features is not None:
        args.features = list(set(args.features) - set(args.ids).union({args.Class, args.eval, args.set}))
        # remove columns that are not feature, id, eval or class
        infile = infile[[k for k in args.features + args.ids + [args.Class, args.eval, args.set] if k in infile]]
        log.info("{}/{} given features found in infile.".format(np.sum(np.in1d(args.features, infile.columns)), len(args.features)))

    # train/test or only test?
    if args.model is None:
        # do the feature prep, training and testing
        if args.cv is not None:
            if args.set in infile:
                log.info("Ignoring set column, using CV.")
                infile.drop(columns=args.set, inplace=True)

            tests, clfs = [], []
            for train, test in cv_samples(infile, args.cv):
                test, clf = traintest(train, test, args.Class, args.threshold, args.n_estimators, args.criterion, args.ids + [args.eval])
                tests.append(test)
                clfs.append(clf)
            test = pd.concat(tests)

        elif args.set not in infile:
                log.info(f"CV and set missing, using all {len(infile)} data points for training.")
                test, clf = traintest(infile, pd.DataFrame({c:[] for c in infile.columns}), args.Class, args.threshold, args.n_estimators, args.criterion, args.ids + [args.eval])
                clfs = [clf]

        else:
            sets = set(s for ss in infile[args.set] for s in ss.split(",") if len(s) > 0)
            if sets == {'train', 'test'}:
                train = infile.loc[infile[args.set] == 'train', :].drop(columns=args.set)
                test = infile.loc[infile[args.set] == 'test', :].drop(columns=args.set)
                log.info("Found train and test set in infile with {} and {} datapoints".format(len(train), len(test)))

                test, clf = traintest(train, test, args.Class, args.threshold, args.n_estimators, args.criterion, args.ids + [args.eval])
                clfs = [clf]

            else:
                assert len(sets) > 1, f"Only {len(sets)} sets found in infile"
                log.info("Found {} sets in infile".format(len(sets)))
                infile[args.set] = [[s for s in ss.split(",") if len(s) > 0] for ss in infile[args.set]]
                tests, clfs = [], []
                for s in sets:
                    testidx = np.asarray([s in ss for ss in infile[args.set]])
                    log.info("Leaving out {} with {} datapoints".format(s, testidx.sum()))
                    train = infile.loc[~testidx, :].drop(columns=args.set)
                    test = infile.loc[testidx, :].drop(columns=args.set)
                    test, clf = traintest(train, test, args.Class, args.threshold, args.n_estimators, args.criterion, args.ids + [args.eval])
                    tests.append(test)
                    clfs.append(clf)
                test = pd.concat(tests)

        if args.out:
            for i, clf in enumerate(clfs):
                joblib.dump(clf, args.out + str(i))

        if test is None: return

    else:
        # testing only
        log.info("Read trained model.")
        clf = joblib.load(args.model)
        log.info("Make prediction on the test set.")
        test = testing(clf, infile, args.ids + [args.Class])

    pred_cols = args.ids + [args.Class, args.eval, "pred"]
    pred_cols = np.intersect1d(pred_cols, test.columns)
    test[pred_cols].to_csv(sys.stdout, sep="\t", index=False)

    if args.importance is not None:
        features = np.setdiff1d(test.columns, args.ids + [args.Class, args.eval, "pred"])
        importance = np.vstack([clf.feature_importances_ for clf in clfs]).mean(axis=0)
        stds = np.std([t.feature_importances_ for clf in clfs for t in clf.estimators_], axis=0)
        importance_table = pd.DataFrame({"feature": features, "importance": importance, "std": stds})
        importance_table.to_csv(args.importance, sep='\t', index=False)

    print_performance(test, args.Class, args.eval, "pred", args.eval_threshold)



if __name__ == '__main__':
    args = get_parser().parse_args()
    args.infile = sys.stdin
    main(args)

debug = False
if debug:
    args = get_parser().parse_args('-i enzyme acceptor source cid -c reaction -e rate'.split(' '))
    args.infile = '/Users/degnbol/biosustain/gt/randomforest/rdkit-desc_muscle-ntern-hmm_median_gtpred_selected/traintest.tsv'

