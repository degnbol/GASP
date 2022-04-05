#!/usr/bin/env zsh
for file in pred*-yield.tsv; do
    echo -n "${file:r}\t"
    mlr --tsv --from $file put '$class = $Yield < 50 ? 0 : 1' |
        rocauc.py -H -p pred -c class /dev/stdin -o ${file:r}.png
done > auc.tsv

# we see that AUC is bad for chem files (novel acceptors)
# but good on seqs files (novel sequences)
# but it seems only good on differentiating yield above below 50%.
# Most have yield=0% and towards the lower end, 
# so it is good at finding very reactive acceptors for new sequences but 
# not very good at differentiating between low and mild reactivity.
