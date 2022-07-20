#!/usr/bin/env zsh
ROOT=`git root`
LIB=$ROOT/tools/degnlib/subfunctions

for file in pred*.gz; do
    gzcat $file | mlr --tsv cut -f cid,enzyme,pred > $file:r
done

grep '^seq' $ROOT/results/*features/selection/species_select.txt > speciesSeqSelect.tmp
INFILES="pred.xlsx.info speciesSeqSelect.tmp"
INFILES="${INFILES} `ls pred*.tsv`"
SHEETS="Overview species_selected_seq_features"
SHEETS="${SHEETS} `for file in pred*.tsv; do echo $file:r | sed 's/pred_//' | sed 's/Amb//'; done`"
echo $INFILES
echo $SHEETS
$LIB/table_excel.py ${=INFILES} -s ${=SHEETS} -o "pred.xlsx" -f


