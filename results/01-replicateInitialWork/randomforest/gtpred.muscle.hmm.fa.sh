hmmalign --trim --amino --outformat A2M gtpred.muscle.hmm `git root`/data/GT-Predict/Active_enzymes_protein_sequences.txt | fasta.py delete-lower > gtpred.muscle.hmm.fa
