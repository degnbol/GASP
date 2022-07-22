#!/usr/bin/env zsh
# TRIM(CLEAN(...)) was used in excel on everything.
# Export table as csv with UTF-8 in excel, to preserve alpha and beta.
# Then convert to tsv and remove zero width spaces, BOM and other annoying things:
sed 's/\xEF\xBB\xBF//g' 20201218_UGT_activitydata_from_lit.csv |
    mlr -c2t cat > 20201218_UGT_activitydata_from_lit.tsv

