#!/usr/bin/env zsh
# INSTALL: brew install fasttree
mlr cut -f enzyme,seq dataset1.tsv | sed 1d | sort -u | sed 's/^/>/' | tr '\t' '\n' | fasttree -slow > dataset1.tsv.tree
