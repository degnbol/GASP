#!/usr/bin/env zsh
ROOT=`git root`
LIB=$ROOT/tools/degnlib/subfunctions
$LIB/table_excel.py betanin-pred.tsv betaninCazy-pred.tsv DON-pred.tsv -s 'betanidin vs tested' 'betanidin vs CAZy' 'DON vs CAZy' -o DON+betanidin-pred.xlsx -fw
