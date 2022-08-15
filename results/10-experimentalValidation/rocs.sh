#!/usr/bin/env zsh
mkdir -p rocs
mlr -t --from blosum62Amb-pred.tsv uniq -f enzyme | sed 1d | while read enz; do
    echo $enz; mlr -t --from blosum62Amb-pred.tsv filter '$enzyme == "'$enz'"' > temp.tsv
    rocauc.py temp.tsv -c Yield -p pred -t $enz -o rocs/$enz.png
done
rm temp.tsv
