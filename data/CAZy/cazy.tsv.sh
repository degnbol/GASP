#!/usr/bin/env zsh
grep -v '^protein_name' CAZy_DB_GlycosylTransferases_11-02-2021.csv | cut -f 1,2,5,6 | cat <(echo "name\tfamily\tspecies\tgenbank") - | mlr --tsv filter '$family == "GT1"' then cut -x -f family | sed 's/^ //' > cazy.tsv
