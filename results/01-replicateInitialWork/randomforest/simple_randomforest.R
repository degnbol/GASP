#!/usr/bin/env Rscript
library(data.table)
library(seqinr)
library(randomForest)
library(pROC)
suppressPackageStartupMessages(library(here))

cwd = function(s) paste0(here(), "/", s)

acc.ints.melt = fread("gtpred_features.tsv", drop="seq")
hts.ugt.melt = fread("hts_ugt_features.tsv", drop="seq")

alignment = read.fasta("gtpred.muscle.hmm.fa", seqtype="AA", whole.header=T, set.attributes=F)
alignment = data.table(enzyme=names(alignment), matrix(unlist(alignment), nrow=length(alignment), byrow=T))
align.melt = melt(alignment, "enzyme", variable.name="pos", value.name="code")[order(enzyme)]

blosum62 = fread(cwd("data/NCBI/blosum62.tsv"))
eye = fread(cwd("data/NCBI/match.tsv"))
setnames(blosum62, "*", "-"); setnames(eye, "*", "-")
setnames(blosum62, "_", "code"); setnames(eye, "_", "code")
blosum62[code=="*", code:="-"]; eye[code=="*", code:="-"]
alphabet = blosum62$code

align.melt.blosum = blosum62[align.melt, on="code"]
align.melt.eye = eye[align.melt, on="code"]
align.melt.blosum[,code:=NULL]
align.melt.eye[,code:=NULL]
align.melt.blosum = melt(align.melt.blosum, c("enzyme", "pos"), variable.name="code")[order(enzyme, pos, code)]
align.melt.eye = melt(align.melt.eye, c("enzyme", "pos"), variable.name="code")[order(enzyme, pos, code)]
align.blosum = dcast(align.melt.blosum[,pos_code:=paste(pos, code, sep="_")], enzyme ~ pos_code, value.var="value")
align.eye = dcast(align.melt.eye[,pos_code:=paste(pos, code, sep="_")], enzyme ~ pos_code, value.var="value")

# for(col in names(alignment)[-1]) alignment[,eval(col):=factor(get(col), levels=alphabet)]
acc.ints.melt = alignment[acc.ints.melt, on="enzyme"]
hts.ugt.melt = alignment[hts.ugt.melt, on="enzyme"]
# acc.ints.melt = align.blosum[acc.ints.melt, on="enzyme"]
# hts.ugt.melt = align.blosum[hts.ugt.melt, on="enzyme"]

# space and dash is problematic for randomForest
names(acc.ints.melt) = sub(" ", "_", names(acc.ints.melt))
names(acc.ints.melt) = sub("-", "_", names(acc.ints.melt))
names(hts.ugt.melt) = sub(" ", "_", names(hts.ugt.melt))
names(hts.ugt.melt) = sub("-", "_", names(hts.ugt.melt))

hts.ugt.melt.rate = hts.ugt.melt[, .(enzyme, acceptor, rate)]

acc.ints.melt$class = factor(acc.ints.melt$class)
acc.ints.melt$Family = factor(acc.ints.melt$Family)
hts.ugt.melt$Family = factor(hts.ugt.melt$Family, levels=levels(acc.ints.melt$Family))
acc.ints.melt$Num_OH = as.numeric(acc.ints.melt$Num_OH)
hts.ugt.melt$Num_OH = as.numeric(hts.ugt.melt$Num_OH)
acc.ints.melt[,acceptor:=NULL]
hts.ugt.melt[,acceptor:=NULL]
acc.ints.melt[,enzyme:=NULL]
hts.ugt.melt[,enzyme:=NULL]
hts.ugt.melt[,rate:=NULL]

forest = randomForest(class ~ ., acc.ints.melt)
evaluation = cbind(hts.ugt.melt.rate, pred=predict(forest, hts.ugt.melt, type="prob")[,"1"])
# evaluation = cbind(hts.ugt.melt.rate, pred=predict(forest, hts.ugt.melt))

cor(evaluation$pred, evaluation$rate)
auc(evaluation$rate > median(evaluation$rate), evaluation$pred)

plot(evaluation$pred, evaluation$rate)

