#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(here))

DT = fread(paste0(here(), "/data/reactions/reactions.tsv"))

# correct conjugate bases
conjugates = fread("conjugates.tsv")
for(i_row in 1:nrow(conjugates)) {
    base = conjugates[i_row, "base"][[1]]
    acid = conjugates[i_row, "acid"][[1]]
    DT[cid==base, cid:=acid]
}

fwrite(DT, "reactions.tsv", sep='\t')

