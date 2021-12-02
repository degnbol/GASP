#!/usr/bin/env zsh
mlr --tsv uniq -f genbank cazy.tsv | sed 1d > GT1.genbank
