#!/usr/bin/env Rscript
library(data.table)
library(seqinr)
library(randomForest)
library(pROC)
suppressPackageStartupMessages(library(here))

setwd(paste0(here(), "/data"))

hts_ugt_rates = fread("../Fatemeh_eval/all-experimental.tsv")
setnames(hts_ugt_rates, "V1", "enzyme")

tmh_rates = fread("raw/pTMH/List_Company_ptMH.tsv")
setnames(tmh_rates, c("compound", "UGT"), c("acceptor", "enzyme"))

acc.props = fread("raw/gtpred/acceptor_interaction_data.txt", select=2:22)
acc.ints = as.matrix(fread("raw/gtpred/acceptor_interaction_data.txt", select=c(2,23:75)), "Name")
# convert away from their weird classification system where 0 = not an enzyme-molecule match, 1= match, anything else is unknown
acc.ints[acc.ints > 1] = NA
acc.ints = data.table(acc.ints, keep.rownames="acceptor")
setnames(acc.props, "Name", "acceptor")
acceptor_families = c('Flavonoid', 'Coumarins', 'Cytokinins', 'Cinnamic acid', 'Benzoate', 'Jasmonic acid', 'Gibberellins', 'Auxin', 'Abscisic acid', 'Sinapic acid', 'Other')
acc.props[, Family:=acceptor_families[Family]]


# then between the two projects
setnames(hts_ugt_rates, "BenzAdenine", "BenzylAdenine")
setnames(hts_ugt_rates, "CinAcid", "CinnamicAcid")

hts.ugt_acceptors = names(hts_ugt_rates)[-1]
acc.ints_acceptors = gsub(" ", "", acc.ints$acceptor)
for(acc in hts.ugt_acceptors) {
    # match on the whole phrase
    idx = grep(paste0("^", acc, "$"), acc.ints_acceptors, ignore.case=T)
    if(length(idx) == 1) {
        message("matching acceptor from HTS UGT: ", acc, " to acc. ints.: ", acc.ints$acceptor[[idx]])
        acc.ints[idx, acceptor:=acc]
        acc.props[idx, acceptor:=acc]
    } else {
        # match on partial phrase
        idx = grep(acc, acc.ints_acceptors, ignore.case=T)
        if(length(idx) == 1) {
            message("matching acceptor from HTS UGT: ", acc, " to GT-Predict: ", acc.ints$acceptor[[idx]])
            acc.ints[idx, acceptor:=acc]
            acc.props[idx, acceptor:=acc]
        }
    }
}


fwrite(acc.props, "acceptor_physiochemical.tsv", sep='\t')


acceptors = hts.ugt_acceptors[hts.ugt_acceptors %in% acc.ints$acceptor]
message("matched on ", paste(acceptors, collapse=", "))

# not mapped:
# hts.ugt_acceptors[!hts.ugt_acceptors %in% acc.ints$acceptor]
# acc.ints$acceptor[!acc.ints$acceptor%in%hts.ugt_acceptors]

acc.ints.melt = na.omit(melt(acc.ints, "acceptor", variable.name="enzyme", value.name="class"))
hts.ugt.melt = melt(hts_ugt_rates, "enzyme", variable.name="acceptor", value.name="rate")
acc.ints.melt = acc.props[acc.ints.melt, on="acceptor"]
hts.ugt.melt = acc.props[hts.ugt.melt[acceptor%in%acceptors], on="acceptor"]

# according to a text file in the code folder zero means no group present and value larger than zero means group is present.
acc.ints.melt[is.na(acc.ints.melt)] = 0
hts.ugt.melt[is.na(hts.ugt.melt)] = 0

# add sequence info
alignment = read.fasta("muscle.fa", seqtype="AA", whole.header=T, set.attributes=F, as.string=T)
alignment = data.table(enzyme=names(alignment), seq=unlist(alignment))
acc.ints.melt = alignment[acc.ints.melt, on="enzyme"]
hts.ugt.melt = alignment[hts.ugt.melt, on="enzyme"]

fwrite(acc.ints.melt, "gtpred_features.tsv", sep="\t")
fwrite(hts.ugt.melt, "hts_ugt_features.tsv", sep="\t")

