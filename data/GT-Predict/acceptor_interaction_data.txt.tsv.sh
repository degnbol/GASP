#!/usr/bin/env zsh
cut -f2,23- acceptor_interaction_data.txt |
    sed '1s/Name/acceptor/' | sed '1s/UGT/At_/g' > acceptor_interaction_data.txt.tsv
