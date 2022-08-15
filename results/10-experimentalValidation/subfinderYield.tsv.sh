#!/usr/bin/env zsh
mlr -t --from SubFinder.tsv rename enzyme,Enzymes,metabolite,smiles + join -j Enzymes,smiles -f forSubFinder.tsv + rename 'Prediction score,pred' > subfinderYield.tsv
