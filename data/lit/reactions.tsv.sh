#!/usr/bin/env zsh

# fix microsoft related characters
function unwin() {
    tr "′’ʹ" "'" | tr '‐' '-' | tr -d ' '
}

unwin < raw/sebastian/gt_all_data.csv > tmp1.csv
unwin < raw/sebastian/gt_acceptors_unique_CIDmanual.csv > tmp2.csv

mlr --icsv --otsv --from tmp1.csv join -j sugar_acceptor -f tmp2.csv then \
    rename 'C-GT,enzyme,sugar_acceptor,acceptor,CID_acceptor,cid,diC-GT_activity,C_reaction,O-GT_activity,O_reaction,kobs_1/s,rate' then \
    put 'if ($C_reaction == "yes" || $O_reaction == "yes") {$reaction = 1} else {$reaction = 0}' then \
    cut -f enzyme,reaction,rate,cid,acceptor then put '$source = "sebastian"' then filter '$cid != ""' > reactions.tsv

rm tmp?.csv

