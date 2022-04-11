#!/usr/bin/env zsh
echo "file\tthres\tn0\tn1\tauc" > auc.tsv
for thres in 25 50 75; do
    for file in pred*-yield.tsv; do
        echo -n "${file:r}\t$thres\t"
        mlr --tsv --from $file put '$class = $Yield < '$thres' ? 0 : 1' > class.tsv.tmp
        mlr --tsv --from class.tsv.tmp uniq -c -f class then cut -f count | sed 1d | tr '\n' '\t'
        rocauc.py -H -p pred -c class -o /dev/null class.tsv.tmp
    done
done | tee -a auc.tsv

rm class.tsv.tmp

# we see that AUC is bad for chem files (novel acceptors)
# but good on seqs files (novel sequences)
# but it seems only good on differentiating yield above below 50%.
# Most have yield=0% and towards the lower end, 
# so it is good at finding very reactive acceptors for new sequences but 
# not very good at differentiating between low and mild reactivity.
