#!/usr/bin/env python3
import os

os.chdir(os.path.expanduser("~/biosustain/gt/feature_importance"))


nFeatures = 0
added = []

with open("feature_importance.sh.e31200019") as infile:
    print("nFeatures\tadded\ttryFeature\tAUC")
    for line in infile:
        line = line.strip()
        if line.startswith("Train random forest model with "):
            nFeatures = int(line.split()[5])
        elif line.startswith("AUC = "):
            tryFeature = line.split()[6]
            auc = float(line.split()[2])
            print(f"{nFeatures}\t{','.join(added)}\t{tryFeature}\t{auc}")
        elif line.startswith("Adding "):
            added.append(line.split()[1])

