#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(webchem))
library(optparse)

if(F) {  # test example
    cids = unique(fread("~/biosustain/gt/data/reactions.tsv")$cid)
    args = list(props="CanonicalSMILES")
} else {
    arg_parser = OptionParser("%prog [options] < infile > outfile", option_list=list(
        make_option(c("-i", "--cid"), type="character", help="Column name with CIDs. Default is assuming the file simply has a CID on each line."),
        make_option(c("-H", "--header"), action="store_true", default=F, help="Indicate that infile has a header. Just a shorthand for --cid cid."),
        make_option(c("-u", "--unique"), action="store_true", default=F, help="Make ids unique."),
        make_option(c("-p", "--props"), type="character", help="List props wanted separated by comma. Default is a custom selection. Set as \"all\" to get all available properties.")
    ))
    args = parse_args(arg_parser)
    if(args$header && is.null(args$cid)) args$cid = "cid"
    if(!is.null(args$cid)) {
        cids = fread("cat /dev/stdin")[[args$cid]]
    } else {
        cids = scan(file("stdin"), quiet=T)
    }
    # ensures that cids are valid, strips spaces at ends
    cids = as.integer(cids)
    if(args$unique) cids = unique(cids)
    message(paste(length(cids), "CIDs given."))
    if(!is.null(args$props)) {
        if (args$props == "all") {
            args$props = NULL
        } else {
            args$props = strsplit(args$props, ',')[[1]]
            args$props[args$props == "smiles"] = "CanonicalSMILES"
        }
    } else {
        args$props = c(
            "CanonicalSMILES",
            "MolecularWeight", 
            # "XLogP", contains NA
            "TPSA", "Complexity", "Charge", "HBondDonorCount", 
            "HBondAcceptorCount", "RotatableBondCount", "HeavyAtomCount", 
            "AtomStereoCount", "BondStereoCount")
    }
}

# pc_prop handles multiple cids and integers even though the help info does not indicate that.
out = pc_prop(cids, args$props)
setnames(out, "CID", "cid")
setnames(out, "CanonicalSMILES", "smiles")
fwrite(out, sep="\t")
