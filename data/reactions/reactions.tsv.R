#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(here))

setwd(paste0(here(), "/data"))

hts_ugt.rates = fread("Fatemeh_eval/all-experimental.tsv")
setnames(hts_ugt.rates, "V1", "enzyme")
hts_ugt.rates = melt(hts_ugt.rates, "enzyme", variable.name="acceptor", value.name="rate")
hts_ugt.rates$source = "HTS_UGT"
tmh.rates = fread("pTMH.tsv", drop=c("sequence", "species"))
tmh.rates$source = "pTMH"
rates = rbind(hts_ugt.rates, tmh.rates)
# based on work in visualize folder and discussion with David we decided to set the threshold for reactivity per enzyme 
# by fitting a normal distribution and looking for outliers.
rates[,pNormEnz:=pnorm(rate, sd=mean(rate), lower.tail=F), by=enzyme]
rates[,qNormEnz:=p.adjust(pNormEnz), by=enzyme]
rates[(pNormEnz < 0.05) & (qNormEnz < 0.05), NormEnz:=1.0]
rates[(pNormEnz > 0.05) & (qNormEnz > 0.05), NormEnz:=0.0]
rates[is.na(NormEnz), NormEnz:=0.5]
# sanity check: how many of each?
rates[,.N,by=NormEnz]
# rm temp
rates[,pNormEnz:=NULL]
rates[,qNormEnz:=NULL]
setnames(rates, "NormEnz", "reaction")

gtpred.reactions = fread("reactions/gtpred_reactions.tsv")
gtpred.reactions$source = "GT-Predict"
DT = merge(gtpred.reactions, rates, all=T, by=c("acceptor", "enzyme", "source", "reaction"))

acc.raw2cid_title = fread("reactions/raw-acceptor_cid_title.tsv", sep='\t')
DT = merge(DT, acc.raw2cid_title, by.x="acceptor", by.y="raw")
DT[, acceptor:=NULL]
setnames(DT, "title", "acceptor")

sebastian = fread("lit/reactions.tsv")
sebastian$source = "sebastian"
sebastian[!is.na(rate), reaction:=1]  # we define lit data where a rate is reported as reactive per discussion with David
DT = rbind(DT, sebastian)
DT[,cid:=as.integer(cid)] # removes spaces at ends and validates

# NOTE that the CIDs are not validated yet, i.e. may be for chemicals that contain CO-. See results/*validat*
fwrite(DT, "reactions/reactions.tsv", sep='\t')

