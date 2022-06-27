#!/usr/bin/env Rscript
library(ggplot2)
library(pROC)
library(data.table)
suppressPackageStartupMessages(library(here))
# go to git root
setwd(here())

# gt.pred = fread("data/GT-Predict/enzyme_interaction_data.txt")
gt.pred = fread("data/reactions/gtpred_reactions_enz.tsv")
acc.int = fread("data/GT-Predict/acceptor_interaction_data.txt", drop=c(1,3:22))
acc.int$Name = sub("^ ", "", acc.int$Name)
acc.int$Name = sub("(?)", "()", acc.int$Name, fixed=T)
acc.int$Name = sub("cinamic acid", "cinnamic acid", acc.int$Name, fixed=T)


length(colnames(gt.pred)[-1])
length(colnames(acc.int)[-1])
length(intersect(colnames(gt.pred), colnames(acc.int)))
length(intersect(gt.pred$acceptor, acc.int$Name))
(unmapped.acceptors = gt.pred[!acceptor%in%acc.int$Name, acceptor])
(unmapped.Names = acc.int[!Name%in%gt.pred$acceptor, Name])


gt.pred.melt = data.table::melt(gt.pred, "acceptor", variable.name="enzyme", value.name="class")
acc.int.melt = data.table::melt(acc.int, "Name", variable.name="enzyme", value.name="acc_class")

DT = merge(gt.pred.melt, acc.int.melt, by="enzyme", allow.cartesian=T)
DT = na.omit(DT)
identities = DT[, .(identities=mean(class==acc_class)), by=c("acceptor", "Name")]
identities[Name%in%unmapped.Names & acceptor%in%unmapped.acceptors & (identities == 1)]

# trans-Zeatin vs Zeatin have identical values in this data. same for 2-Methyl-umbelliferone vs 4-Methyl-umbelliferone
# searching for zeatin gives tans-zeatin as first result. It might be more common that the cis version
# searching for 4-Methyl-umbelliferone gives many results as opposed to 2-Methyl-umbelliferone
# looked on https://pubchem.ncbi.nlm.nih.gov/#query=zeatin and https://www.ebi.ac.uk/chebi/advancedSearchFT.do

# all enzyme from gt.pred are found in acc.int so gt.pred is just a subset of the data in acc.int
names(gt.pred)[!names(gt.pred)%in%names(acc.int)]


