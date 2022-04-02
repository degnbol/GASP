#!/usr/bin/env zsh
cut -f2,23- raw/gtpred/acceptor_interaction_data.txt |
    sed '1s/Name/acceptor/' | sed '1s/UGT/At_/g' > gtpred_reactions_acc.tsv
