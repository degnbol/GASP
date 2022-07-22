#!/usr/bin/env zsh
# remove two letter species indication, then remove accession so the ones in the manual file is used instead.
mlr --tsv --from 20201218_UGT_activitydata_from_lit.tsv put '$Protein = substr($Protein, 2, -1)' then \
    cut -x -f Accession then \
    join -j Protein -f 20201218_UGT_activitydata_from_lit.manual.tsv > 20201218_UGT_activitydata_from_lit.AA.tsv
