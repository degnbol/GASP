#!/usr/bin/env zsh
mlr join -j enzyme,cid -f matchAmb-pred.tsv --lk pred --lp match_ + \
    rename pred,blosum62_pred blosum62Amb-pred.tsv |
    auc -p match_pred blosum62_pred -c Yield -t 50 -g enzyme -n |
    mlr filter '$Yield_thres_50_1 > 4' + sort -nr Yield_thres_50_blosum62_pred_AUC > enzymes.auc.tsv
