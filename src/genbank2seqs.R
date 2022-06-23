#!/usr/bin/env Rscript
# only implemented for protein currently
# USAGE: genbank2seqs.R < ids.genbank > seqs.faa
library(rentrez)
ids = readLines("stdin")
n_ids = length(ids)
message("#ids = ", n_ids)
chunk_size = 200
for(i in seq(1, n_ids, chunk_size)) {
    message(sprintf("%.2f", i/n_ids*100), "%")
    cat(entrez_fetch("protein", ids[i:min(i+chunk_size-1, n_ids)], rettype="fasta", email="christian.degnbol@gmail.com", api_key="1f55b5052ac1a0ef0fd8869f1a945d6bcd09"))
}
