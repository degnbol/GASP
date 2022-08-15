#!/usr/bin/env zsh
auc pred_all.tsv -c Yield -t 50 -g 'Polyphenol Name' -p pred -n -b -a -C > pred_all-aucChem.tsv
# add corr
mlr -t --from pred_all.tsv stats2 -a corr -f Yield,pred -g 'Polyphenol Name' + join -j 'Polyphenol Name' -f pred_all-aucChem.tsv > temp && mv temp pred_all-aucChem.tsv
