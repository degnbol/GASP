#!/usr/bin/env zsh
# curate most from EMBL
`git root`/src/acc2aa.sh UGT_reactivity_Natalia.tsv
# some were looked up on uniprot and NCBI manually.
# Do a left join where we annotate the manual ones and
# emit unpaired from left file --ul to get the rest.
mlr --tsv --from UGT_reactivity_Natalia.AA.manual.tsv cut -x -f Source then \
    join --ul -j Accession -f UGT_reactivity_Natalia.AA.tsv > temp && mv temp UGT_reactivity_Natalia.AA.tsv
