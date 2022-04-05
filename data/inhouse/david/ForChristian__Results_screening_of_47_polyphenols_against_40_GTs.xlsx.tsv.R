#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

bname = "ForChristian__Results_screening_of_47_polyphenols_against_40_GTs.xlsx"

# read
DT = fread(paste0(bname, ".yields.tsv"), header=TRUE)
DT.enz = fread(paste0(bname, ".enzymes.tsv"))
DT.cid = fread("20220215__forML.xlsx.cas.tsv")

DT.melt = melt(DT, id.vars=c("Polyphenol Name", "CAS Number"), variable.name="ID", value.name="Yield", variable.factor=FALSE, value.factor=FALSE)
# made a factor, we convert so they are the same type.
DT.melt[,ID:=as.integer(ID)]
DT.join = DT.enz[DT.melt, on="ID"]

# annotate CIDs
setnames(DT.join, "CAS Number", "cas")
DT.join.cid = DT.cid[DT.join, on="cas"]

fwrite(DT.join.cid, paste0(bname, ".tsv"), sep='\t', quote=FALSE)
