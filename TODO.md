[x] are we following DOME? see Ditte email.
 - we are partly, but if we want to address and cover everything mentioned in 
   the DOME paper then I have a lot more work to do.
x.] look thru writing, are they saying something silly about ML?
[x] make phylo of dataset 1, 1031 datapoints maybe with clustal omega (what 
david used for other phylo). Show chems on other axis. show bool reactivity 
with red, gray, green and put raw rates in fig as numbers.
[x] dataset as a function of phylo
[.] highest prio. compare performance vs others: GT-predict, etc. in results/16-cmpPerform.
  [x] ESP. Worse than random performance, even on metabolites in its trainset.
  [x] GT-Predict can be next. It isn't pan-specific so we can only predict on 
      chems in GT-Predict that also happens to be in dataset 1 and then any 
      sequence in dataset 1. RESULT: great performance, obvious on identical 
      seqs, but almost .7 on seqs with neighbor.
  [ ] pan-specific randomforest. Has to do better than .7.
[ ] lit data moves to eval.
[x] neg feature selection. - Running on spartan.
    Takes a lot of time for seq select so only decision tree for that.
    Not really improving performance.
