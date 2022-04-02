#!/usr/bin/env zsh
grep Add feature_importance.sh.e33141314 | cut -d' ' -f2- | sed 's/, /\n/' | sort -u > speciesTop_select.txt
