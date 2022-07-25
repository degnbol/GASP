#!/usr/bin/env zsh
cd `git root`/data/reactions
./gtpred_reactions_acc.tsv.sh
./gtpred_reactions_enz.tsv.sh
./gtpred_reactions.tsv.R
./rawAcceptor2cid.tsv.R
./rawAcceptor_cid_title.tsv.sh
./reactions.tsv.R
