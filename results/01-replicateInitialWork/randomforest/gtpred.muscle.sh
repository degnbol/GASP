cat `git root`/GT-Predict/Active_enzymes_protein_sequences.txt | muscle -quiet 2> gtpred.muscle.log | fasta.py unwrap > gtpred.muscle.fa
