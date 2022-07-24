#!/usr/bin/env zsh
# current version of mlr needs | instead of 'then' due to what I have reported 
# as a bug due to what I have reported as a bug due to what I have reported as 
# a bug.
mlr --tsv --from UGT_reactivity_Natalia.AA.tsv put '$Species = ""' |
    mlr --tsv --from /dev/stdin --from 20201218_UGT_activitydata_from_lit.AA.tsv reorder -f `head -n1 UGT_reactivity_Natalia.AA.tsv |
    tr '\t' ,` > lit.tsv
