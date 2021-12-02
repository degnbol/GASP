#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(webchem))
library(optparse)

if(F) {  # test example
    SMILESs = "CC"
} else {
    arg_parser = OptionParser(
        description="Convert SMILES to CID. Output table has columns smiles,cid.",
        "%prog [options] < infile > outfile", option_list=list(
        make_option(c("-i", "--smiles"), type="character", help="Column name with SMILESs. Default is assuming the file simply has a SMILES on each line."),
        make_option(c("-H", "--header"), action="store_true", default=F, help="Indicate that infile has a header. Just a shorthand for --smiles smiles."),
        make_option(c("-u", "--unique"), action="store_true", default=F, help="Make ids unique.")
    ))
    args = parse_args(arg_parser)
    if(args$header && is.null(args$smiles)) args$smiles = "smiles"
    if(!is.null(args$smiles)) {
        SMILESs = fread("cat /dev/stdin")[[args$smiles]]
    } else {
        SMILESs = fread("cat /dev/stdin", header=F)[["V1"]]
    }
    if(args$unique) SMILESs = unique(SMILESs)
}

out = get_cid(SMILESs, from="smiles", match="first")
setnames(out, "query", "smiles")
fwrite(out, sep="\t")
