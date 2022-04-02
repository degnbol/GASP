#!/usr/bin/env zsh
LC_ALL=C sed $'s/[^[:print:]\t]//g' raw/gtpred/enzyme_interaction_data.txt |
    cut -f2- | sed 2d | tr -d '"' | sed '1s/^/acceptor/' | sed '1s/UGT/At_/g' > gtpred_reactions_enz.tsv
