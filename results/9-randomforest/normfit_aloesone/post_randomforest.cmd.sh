#!/usr/bin/env zsh
# after running randomforest training on GT-Predict + in-house + lit and prediction on aloesone etc.

# cleanup meaningless cols
mlr -I --tsv cut -x -f "rate,reaction,source" pred*.tsv

for file in pred*.tsv; do
    # make top100
    mlr --tsv --from $file sort -nr pred | mlr --tsv head -n 100 -g acceptor | mlr --tsv sort -f acceptor > top100_${file:r}.tsv
    # annotate data source
    mlr --tsv --from `git root`/results/*generateNegatives/reactions.tsv uniq -f "enzyme,source" then join -j enzyme -f $file then cut -x -f rate,reaction > source_${file:r}.tsv
    # combine into excel
    table.py excel top100_${file:r}.tsv source_${file:r}.tsv -s "cazy top 100" "with source" -fw -o "aloesoneEtc_${file:r:t}.xlsx"
done

