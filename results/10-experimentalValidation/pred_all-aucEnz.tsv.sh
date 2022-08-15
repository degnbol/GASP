#!/usr/bin/env zsh
auc pred_all.tsv -c Yield -t 50 -g enzyme -p pred -n -b -a -C > pred_all-aucEnz.tsv
# add corr
mlr -t --from pred_all.tsv stats2 -a corr -f Yield,pred -g enzyme + join -j enzyme -f pred_all-aucEnz.tsv > temp && mv temp pred_all-aucEnz.tsv
