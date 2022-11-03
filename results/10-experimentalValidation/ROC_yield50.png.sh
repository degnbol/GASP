#!/usr/bin/env zsh
mlr -t --from smpl_blosum62Amb-pred.tsv put '$yield50 = int($Yield > 50)' > temp.tsv
rocauc.py -c yield50 -p pred temp.tsv -o ROC_yield50.png
rm temp.tsv
