#!/usr/bin/env zsh
./positives.tsv.sh
./gen_negatives.py
./negatives.tsv.sh
`git root`/data/reactions-neg.tsv.sh
