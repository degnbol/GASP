#!/usr/bin/env zsh
ROOT=`git root`

# it changes delim from comma to semicolon, 
# renames the output columns and 
# removes our mapping columns so we rename them back to join tables.
mlr -c --ifs \; cat + rename metabolite,Metabolites,enzyme,Enzymes 5*_output.csv > ESP_output.csv

# assert that all inputs were valid
if [[ `mlr -c --ho --from ESP_output.csv uniq -f 'valid input' | xargs` != "True" ]]; then
    echo "Not all inputs were valid."
    exit 1
fi

reactions=`\ls $ROOT/results/*generateNegatives/reactions.tsv`

mlr --c2t join -f ESP_output.csv -j Metabolites,Enzymes +\
    rename '#metabolite in training set,nMetabInESP,Prediction score,pred_ESP' +\
    cut -f nMetabInESP,pred_ESP,cid,enzyme ESP_input_*.csv |
    mlr -t join -f $reactions -j cid,enzyme --lk reaction + uniq -a > ESP_output.tsv

rm ESP_output.csv

