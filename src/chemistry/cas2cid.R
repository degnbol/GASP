#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(webchem))

CASs = readLines("stdin")

queries = cts_convert(CASs, from="CAS", to="pubchem CID", match="first")

cat("cas\tcid\n")
for(CAS in names(queries)) {
    CID = queries[[CAS]]
    if(is.na(CID)) {
        # fall back on general match since CAS is often listed as a synonym
        CID = get_cid(CAS, match="first")[["cid"]]
    }
    cat(paste0(CAS, '\t', CID, '\n'))
}


