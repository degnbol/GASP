#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(webchem))
suppressPackageStartupMessages(library(here))
library(optparse)


if(F) {  # test example
    cids = unique(fread(paste0(here(), "/data/reactions/reactions.tsv"))$cid)
    args = list(props="CanonicalSMILES")
} else {
    arg_parser = OptionParser(
        usage="%prog [options] < INFILE with CIDs > OUTFILE.tsv", 
        description="Convert CIDs to SMILES using pubchem. Additional chemical properties can also be fetched. Any unrecognized inputs will be assumed to be SMILES and will passthrough to the smiles column.",
        option_list=list(
            make_option(c("-i", "--cid"), type="character", help="Column name with CIDs. Default is assuming the file simply has a CID on each line."),
            make_option(c("-H", "--header"), action="store_true", default=F, help="Indicate that infile has a header. Just a shorthand for --cid cid."),
            make_option(c("-u", "--unique"), action="store_true", default=F, help="Make ids unique."),
            make_option(c("-p", "--props"), type="character", help="List props wanted separated by comma. Default is a custom selection. Set as \"all\" to get all available properties.")
        )
    )
    args = parse_args(arg_parser)
    if(args$header && is.null(args$cid)) args$cid = "cid"
    if(!is.null(args$cid)) {
        cids = fread("cat /dev/stdin")[[args$cid]]
    } else {
        cids = fread("cat /dev/stdin", header=FALSE)[["V1"]]
    }
    if(args$unique) cids = unique(cids)
    message(paste(length(cids), "CIDs/SMILESs given."))
    if(!is.null(args$props)) {
        if (args$props == "all") {
            args$props = NULL
        } else {
            args$props = strsplit(args$props, ',')[[1]]
            args$props[args$props == "smiles"] = "CanonicalSMILES"
        }
    } else {
        args$props = "CanonicalSMILES"
    }
}

# assume cid can be parsed to int and SMILES can't.
# Use strtoi instead of as.integer to avoid warning when unable to parse.
isSMILES = is.na(strtoi(cids))
# cids named UNKNOWN# are assigned since cids are used as IDs through the chemical pipeline to distinguish chemicals.
# Note that some of the downstream code doesn't allow even an underscore in the cid, but UNKNOWN# is tested to work.
out = data.table(cid=paste0("UNKNOWN", 1:sum(isSMILES)), smiles=cids[isSMILES])
if(sum(!isSMILES) > 0){
    # pc_prop handles multiple CIDs and integers even though the help info does not indicate that.
    pub = data.table(pc_prop(cids[!isSMILES], args$props))
    setnames(pub, "CID", "cid")
    setnames(pub, "CanonicalSMILES", "smiles")
    # NOTE: SMILES input entries are placed at the top, so row order is not maintained for mixed input.
    out = rbind(out, pub)
}
# NA is written as empty cell by default.
fwrite(out, sep="\t")

