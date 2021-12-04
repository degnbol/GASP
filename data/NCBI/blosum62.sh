#!/usr/bin/env zsh
# version with and without B Z X since they are not the standard 20 amino acids
egrep -v '^#|^B|^Z|^X' raw/BLOSUM62 | sed '1s/^ /_/' | mlr --ipprint --otsv cut -x -f B,Z,X > blosum62.tsv
egrep -v '^#|^B|^Z|^X' raw/BLOSUM62 | sed '1s/^ /_/' | mlr --ipprint --otsv cat > blosum62Amb.tsv

