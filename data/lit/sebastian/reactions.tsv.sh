#!/usr/bin/env zsh
mlr --icsv --otsv --from `git root`/data/raw/lit/sebastian/gt_all_data.csv join -j sugar_acceptor -f `git root`/data/raw/lit/sebastian/gt_acceptors_unique_CIDmanual.csv then \
    rename 'C-GT,enzyme,sugar_acceptor,acceptor,CID_acceptor,cid,diC-GT_activity,C_reaction,O-GT_activity,O_reaction,kobs_1/s,rate' then \
    put 'if ($C_reaction == "yes" || $O_reaction == "yes") {$reaction = 1} else {$reaction = 0}' then \
    cut -f enzyme,reaction,rate,cid,acceptor then put '$source = "sebastian"' then filter '$cid != ""' > reactions.tsv
