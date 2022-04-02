#!/usr/bin/env zsh
grep Add feature_importance.sh.e33141314 feature_importance.sh.e33138379 feature_importance.sh.e33100432 feature_importance.sh.e33054179 | cut -d' ' -f2- | sed 's/, /\n/' | sort -u > species_select.txt
