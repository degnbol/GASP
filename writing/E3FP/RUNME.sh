#!/usr/bin/env zsh
`git root`/src/chemistry/e3fp-features.py `git root`/results/05-chemicalFeatures/acceptors-props/E3FP.fpz -ck 3 -d E3FP.mat > MDS3.tsv
./plotting.py
