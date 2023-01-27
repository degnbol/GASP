#!/usr/bin/env zsh
mlr -t --from ../../*-dataset1/dataset1.tsv filter '$reaction != 0.5' > testset.tsv
