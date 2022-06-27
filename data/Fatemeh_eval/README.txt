The dataset HTS_UGT contain reaction measurements for a number of enzymes 
and acceptors collected in all-experimental.tsv. All enzymes from this set were given to the 
enzyme specific predictor of GT-Predict, which found nearest neighbor enzyme 
in their dataset and provided the reactivity booleans for all acceptors in 
GT-Predict. These booleans are collected in all-gtpredict.csv.
Then evaluate-gtpredict.R can be run which compares the values for the 
acceptors these two datasets have in common. The performance is quite good 
but this is expected since it's not really predicting on new unseen data, 
rather it is just copying values over for the same already measured 
acceptors, measured on slightly different enzymes. This result mainly shows 
that the GT enzymes will have comparable binding preference if they are 
close in sequence.
