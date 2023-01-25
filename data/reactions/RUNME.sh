#!/usr/bin/env zsh
cd `git root`/data/reactions
./rawAcceptor2cid.tsv.R
./rawAcceptor_cid_title.tsv.sh
./reactions.tsv.R
./reaction-conflicts.tsv.sh
