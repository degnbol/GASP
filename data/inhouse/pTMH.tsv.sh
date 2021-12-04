#!/usr/bin/env zsh
mlr --tsv --from raw/pTMH/List_Company_ptMH.tsv put '$plasmid = "pTMH" . $UGT' then join -j plasmid -f raw/pTMH/Protein_info_pTMH_extract.tsv then \
    rename compound,acceptor then cut -x -f UGT,plasmid > pTMH.tsv
