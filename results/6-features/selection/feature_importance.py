#!/usr/bin/env python3
import argparse
import sys
import numpy as np
import pandas as pd
from multiprocess import Pool
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score as AUC
import logging as log
log.basicConfig(format="%(message)s", level=log.INFO, force=True)


def get_parser():
    parser = argparse.ArgumentParser(
        description="Train a Random Forest classifier with a train and test set and make prediction on the test set. "
                    "Optionally use CV to split infile into train and test. "
                    "Reads a tab-separated table header from stdin and writes prediction table to stdout. "
                    "Columns: id, features and class. Class is 1 for a positive and 0 for a negative.")

    parser.add_argument("-i", "--id", dest="ids", nargs="+", help="Name of columns that are neither feature nor class so should be identifier columns. They are simply copied to output for identification.")
    parser.add_argument("-t", "--identity", dest="threshold", nargs="+", type=float, default=1, help="Threshold for train vs test identity. Can be multiple values to have different thresholds for different groups of features (use -F/--features).")
    parser.add_argument("-o", "--out",
                        help="The model can be saved to this file (e.g. randomforest.joblib). By default the model is not saved.")
    parser.add_argument("-n", "--estimators", dest="n_estimators", type=int, default=10, help="Number of trees in the forest.")
    parser.add_argument("--crit", dest="criterion", default="gini", help="Randomforest criterion for split quality.", choices=["gini", "entropy"])
    # uppercase Class since class is a reserved word
    parser.add_argument("-c", "--class", dest="Class", help="Name of column with class labels in training data. Default=class.", default="class")
    parser.add_argument("-e", "--eval", default="class",
                        help="Name of column with evaluation value in test data. Can be class like infile but floats are also supported. Default=class.")
    parser.add_argument("-F", "--features", nargs="+", help="Define which column names are features. Default is all that are not given as ids, class or eval. "
                                                            "Separate with comma or spaces. If using multiple -t/--identity then spaces separate groups and commas separate features within groups.")
    parser.add_argument("-S", "--select", nargs="+", help="Preselect these features. For -r/--rev this arg selects the features that should always be used, i.e. they will be protected from removal.")
    parser.add_argument("-l", "--log", help="Log file.")
    parser.add_argument("-s", "--set", help="Name of column with set indication. Can have comma separated values. The test set is each unique name found. Default=\"set\".", default="set")
    parser.add_argument("--nproc", type=int, help="Number of cores for multi processing.")
    parser.add_argument("-r", "--rev", action="store_true", help="Remove features one by one instead of adding them.")

    return parser


def get_args(args_string=None):
    args = get_parser().parse_args(None if args_string is None else args_string.split())
    args.ids = [i for i in args.ids if i not in [args.set, args.Class, args.eval]]
    if args.features is not None: args.features = [fs.split(',') for fs in args.features]
    return args

# functions


def top_pred_accuracy(clas, pred):
    """
    As an alternative to AUC, where the focus is on performance for the top scorers.
    :param clas: boolean or int vector with 1s and 0s. True is the class we want to predict with a high prediction score.
    :param pred: float vector. Same length as clas.
    :return float scalar [0;1] or nan if there are no positives
    """
    assert len(clas) == len(pred)
    assert pred.dtype == float
    # avoid pandas
    clas = np.asarray(clas, dtype=bool)
    pred = np.asarray(pred, dtype=float)
    n_pos = clas.sum()
    if n_pos == 0: return np.nan
    # what is the class for the top predicted n_pos datapoints? In order from highest pred score.
    top_pred_class = clas[np.argsort(-pred)[:n_pos]]
    # allocate scoring available
    scoring = np.arange(n_pos, 0, -1)
    # sum from scoring if it is a hit. Normalize so the maximum possible score is 1.
    return scoring[top_pred_class].sum() / scoring.sum()
    

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


def identity_thresholds_idx(train, test, thresholds):
    """
    Get index of rows in train that satisfies a threshold of allowed identity with rows from test.
    :param train: pandas with a data point per row
    :param test: pandas with a data point per row
    :param thresholds: dict mapping from each column name to floats [0;1]
    :return: logical vector
    """
    thresholds = np.asarray([thresholds[c] for c in train])
    train_mat = np.asarray(train)
    test_mat = np.asarray(test)
    return np.asarray([np.mean((row == test_mat) / thresholds, axis=1).max() for row in train_mat]) <= 1


def training(train, Class, n_estimators, criterion):
    """
    Perform training of random forest classifier with proper features given.
    :param train: training features.
    :param Class: name of column with class label.
    :param n_estimators:
    :param criterion:
    :return: trained RandomForestClassifier
    """
    clf = RandomForestClassifier(n_estimators=n_estimators, n_jobs=n_estimators, criterion=criterion)
    clf.fit(train.drop(columns=Class), train[Class])
    return clf


def testing(clf, test, ignore):
    return clf.predict_proba(test.drop(columns=ignore, errors='ignore'))[:, clf.classes_ == 1].squeeze()


def traintest(train, test, Class, n_estimators, criterion, ignore=None, selected=None, threshold=None, nproc=None, rev=False):
    """

    :param train:
    :param test:
    :param Class:
    :param n_estimators:
    :param criterion:
    :param ignore:
    :param selected: Preselected features that must be given to the random forest, i.e. preselection for addition mode and protection for reverse (removal) mode.
    :param threshold: threshold for removing training data that have too high identity with a test data point.
    :param nproc: number of processors. Default: auto-detect
    :param rev: remove features one-by-one instead of adding them.
    :return:
    """
    if ignore is None: ignore = []
    elif type(ignore) is set: ignore = list(ignore)
    elif type(ignore) is not list: ignore = [ignore]
    multi_set = test is None  # then train should be a list of datasets that each can be test set and the remaining trainset
    if multi_set:
        datasets = [ds.drop(columns=ignore, errors='ignore') for ds in train]
        feature_names = set(datasets[0].columns) - {Class}
    else:
        train = train.drop(columns=ignore, errors='ignore')
        feature_names = set(train.columns) - {Class}

    if selected is None: selected = set()
    if rev:
        to_consider = np.asarray(list(feature_names - set(selected)) # consider these for removal
        selected = feature_names

    if type(threshold) == dict:
        _identity_thres = lambda _train, _test: identity_thresholds_idx(_train, _test, threshold)
    elif threshold == 1:
        _identity_thres = None
    else:
        _identity_thres = lambda _train, _test: identity_threshold_idx(_train, _test, threshold)

    # metrics of performance

    def get_performance(train, test, feature_names):
        _train = train[feature_names + {Class}]
        _test = test[feature_names]
        if _identity_thres is not None: _train = _train[_identity_thres(train[feature_names], _test)]
        if (len(train)-len(_train))/len(train) > 0.1:
            log.info(f'Not training on {len(train)-len(_train)}/{len(train)} ({(len(train)-len(_train))/len(train)*100:.2f}%) of training data points due to high identity.')
            if (len(train) - len(_train)) / len(train) > 0.5:
                log.warning(f'Discarded majority. Not training.')
                return np.nan
        clf = training(_train, Class, n_estimators, criterion)
        pred = testing(clf, _test, ignore + [Class])
        # return AUC(test[Class], pred) # instead of AUC we are trying a custom metric biased for top predicted.
        return top_pred_accuracy(test[Class], pred)

    def get_performance_multi(datasets, feature_names):
        performances = []
        for choice in range(len(datasets)):
            train = pd.concat([datasets[i] for i in range(len(datasets)) if i != choice])
            performances.append(get_performance(train, datasets[choice], feature_names))
        return np.nanmean(performances)

    if rev:
        def _get_score(feature_name):
            if multi_set: return feature_name, get_performance_multi(datasets, selected - {feature_name})
            else: return feature_name, get_performance(train, test, selected - {feature_name})
    else:
        def _get_score(feature_name):
            if multi_set: return feature_name, get_performance_multi(datasets, selected + {feature_name})
            else: return feature_name, get_performance(train, test, selected + {feature_name})

    def _get_two_score(feature_name_0):
        unselected_1 = unselected - {feature_name_0}
        np.random.shuffle(unselected_1)
        # no result return values
        score, feature_name_1 = 0, None
        for attempts, feature_name_1 in enumerate(unselected_1):
            if attempts >= 3: break
            if multi_set:
                score = get_performance_multi(datasets, selected + {feature_name_0, feature_name_1})
            else: score = get_performance(train, test, selected + {feature_name_0, feature_name_1})
        return feature_name_0, feature_name_1, score

    def _consider_add():
        consider = set()
        # train for each unselected feature
        unselected = np.asarray(list(feature_names - selected))
        np.random.shuffle(unselected)
        for feature_name, score in pool.imap_unordered(_get_score, unselected):
            if score >= best_last:
                log.info(f"score = {score:.3f} when adding feature {feature_name}")
                if score > best:
                    best = score
                    consider = {feature_name}

        if len(consider) == 0:
            log.info("No single feature improved performance when added. Trying two features.")
            np.random.shuffle(unselected)
            for feature_name_0, feature_name_1, score in pool.imap_unordered(_get_two_score, unselected):
                if score >= best_last:
                    log.info(f"performance = {score:.3f} when adding features {', '.join([feature_name_0, feature_name_1])}")
                    if score > best:
                        best = score
                        consider = {feature_name_0, feature_name_1}

        if len(consider) == 0:
            log.info("Nothing was an improvement. Trying again.")
            return
        log.info("Adding " + ', '.join(consider))
        selected += consider
        best_last = best

    def _consider_remove():
        consider = {}
        # train for each to_consider feature
        np.random.shuffle(to_consider)
        for feature_name, score in pool.imap_unordered(_get_score, to_consider):
            if score >= best_last:
                log.info(f"score = {score:.3f} when removing feature {feature_name}")
                if score > best:
                    best = score
                    consider = {feature_name}

        if len(consider) == 0:
            log.info("No single feature improved performance when removed. Trying two features.")
            np.random.shuffle(to_consider)
            for feature_name_0, feature_name_1, score in pool.imap_unordered(_get_two_score, to_consider):
                if score >= best_last:
                    log.info(f"performance = {score:.3f} when removing features {', '.join([feature_name_0, feature_name_1])}")
                    if score > best:
                        best = score
                        consider = {feature_name_0, feature_name_1}

        if len(consider) == 0:
            log.info("Nothing was an improvement. Trying again.")
            return
        log.info("Adding " + ', '.join(consider))
        selected consider
        best_last = best

    _consider = _consider_remove if rev else _consider_add

    best = best_last = 0.5
    with Pool(nproc) as pool:
        log.info(f"num proc = {pool._processes}")
        for n_features in range(len(selected)+1, len(feature_names)+1)[::-1 if rev else 1]:
            log.info(f"Train random forest model with {n_features:d} feature(s).")
            _consider()

    return selected



def main(args):
    if args.log is not None: log.basicConfig(handlers=[log.StreamHandler(), log.FileHandler(args.log)])

    # reading
    infile = pd.read_table(args.infile)
    if infile.isna().any().any():
        log.warning("NA found in infile in column(s): " + ",".join(infile.columns[infile.isna().any()]))


    # features
    if args.features is None: args.features = [infile.columns]
    args.features = [[f for f in fs if f and f not in args.ids + [args.set, args.Class, args.eval]] for fs in args.features]
    if args.select is not None:
        n_before = len(args.select)
        args.select = [s for s in args.select if s in [f for fs in args.features for f in fs]]
        log.info(f"{len(args.select)}/{n_before} preselected features found.")
    else:
        log.info("No preselected features.")
    if len(args.threshold) == 1:
        args.features = [f for fs in args.features for f in fs]
        log.info(f'Using {len(args.features)} features with threshold {args.threshold}.')
    else:
        assert len(args.threshold) == len(args.features), f"{len(args.threshold)} != {len(args.features)}"
        for t, fs in zip(args.threshold, args.features):
            log.info(f'Using {len(fs)} features with threshold {t}.')

    # make thresholds into dict if there are multiple feature groups
    if len(args.threshold) == 1:
        threshold, = args.threshold
    else:
        threshold = {}
        for t, fs in zip(args.threshold, args.features):
            for f in fs: threshold[f] = t


    if args.set in infile:
        sets = set(s for ss in infile[args.set] for s in ss.split(",") if len(s) > 0)
        if sets == {'train', 'test'}:
            train = infile.loc[infile[args.set] == 'train', :].drop(columns=args.set)
            test = infile.loc[infile[args.set] == 'test', :].drop(columns=args.set)
            log.info(f"Found train and test set in infile with {len(train)} and {len(test)} datapoints")
            print(traintest(train, test, args.Class, args.n_estimators, args.criterion, args.ids + [args.eval], 
                selected=args.select, threshold=threshold, nproc=args.nproc, rev=args.rev))
        else:
            log.info(f"Found {len(sets)} sets in infile")
            infile[args.set] = [[s for s in ss.split(",") if len(s) > 0] for ss in infile[args.set]]
            datasets = []
            for s in sets:
                testidx = np.asarray([s in ss for ss in infile[args.set]])
                log.info(f"Leaving out {s} with {testidx.sum()} datapoints")
                datasets.append(infile.loc[testidx, :].drop(columns=args.set))
            print(traintest(datasets, None, args.Class, args.n_estimators, args.criterion, args.ids + [args.eval],
                selected=args.select, threshold=threshold, nproc=args.nproc, rev=args.rev))
    else:
        raise NotImplementedError("CV and set missing.")


if __name__ == '__main__':
    log.info(" ".join(sys.argv))
    args = get_args()
    args.infile = sys.stdin
    main(args)

debug = False
if debug:
    args = get_args(
        "~/biosustain/gt/features/rdkit-desc_muscle-ntern-hmm_median_notgtpred.tsv --aa seq -c reaction -i enzyme acceptor cid source -t 0.9 -e rate --aa-encoding ~/biosustain/data/ncbi/match.tsv --importance ~/biosustain/gt/permutation_importance/importance.tsv")
    args.infile = "~/biosustain/gt/features/rdkit-desc_muscle-ntern-hmm_median_gtpred.tsv"
    main(args)
