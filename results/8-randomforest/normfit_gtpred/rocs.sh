#!/usr/bin/env zsh
rocauc.py pred_onehot.tsv -p pred -c reaction -t "normfit onehot" -o roc_normfit_onehot.pdf
rocauc.py pred_blosum.tsv -p pred -c reaction -t "normfit blosum62" -o roc_normfit_blosum.pdf
