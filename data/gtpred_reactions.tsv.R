#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(here))

# just a script to decide between enzyme and acceptor file from GT-Predict and make schema better. 
# ids are still raw (human readable at the end of this script)

setwd(paste0(here(), "/data"))

gtpred_enz.reactions = as.matrix(fread("gtpred_reactions_enz.tsv"), "acceptor")
gtpred_acc.reactions = as.matrix(fread("gtpred_reactions_acc.tsv"), "acceptor")
# convert away from their weird classification system where 0/NA = no reaction, 1= reaction, 2= unclear, 3= missing
gtpred_enz.reactions[is.na(gtpred_enz.reactions)] = 0
gtpred_acc.reactions[is.na(gtpred_acc.reactions)] = 0
gtpred_enz.reactions[gtpred_enz.reactions == 2] = 0.5
gtpred_acc.reactions[gtpred_acc.reactions == 2] = 0.5
gtpred_enz.reactions[gtpred_enz.reactions == 3] = NA
gtpred_acc.reactions[gtpred_acc.reactions == 3] = NA
gtpred_enz.reactions = data.table(gtpred_enz.reactions, keep.rownames="acceptor")
gtpred_acc.reactions = data.table(gtpred_acc.reactions, keep.rownames="acceptor")
gtpred_enz.reactions = na.omit(data.table::melt(gtpred_enz.reactions, "acceptor", variable.name="enzyme", value.name="reaction"))
gtpred_acc.reactions = na.omit(data.table::melt(gtpred_acc.reactions, "acceptor", variable.name="enzyme", value.name="reaction"))


enz_acc.cmp = merge(gtpred_enz.reactions, gtpred_acc.reactions, by=c("acceptor", "enzyme"))
cor(enz_acc.cmp$reaction.x, enz_acc.cmp$reaction.y)
# as expected the files are identical when there is overlap in identification, so we use the largest of the two
nrow(gtpred_enz.reactions)
nrow(gtpred_acc.reactions)

gtpred_enz.reactions = gtpred_enz.reactions[order(acceptor, enzyme)]
fwrite(gtpred_enz.reactions, "gtpred_reactions.tsv", sep='\t')
