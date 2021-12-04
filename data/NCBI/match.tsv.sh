#!/usr/bin/env zsh
# without B Z X since they are not the standard 20 amino acids
egrep -v '^#|^B|^Z|^X' raw/MATCH | sed '1s/^ /_/' | mlr --ipprint --otsv cut -x -f B,Z,X > match.tsv
