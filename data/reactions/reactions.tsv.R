#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(here))

setwd(paste0(here(), "/data"))

gtpred.reactions = fread("GT-Predict/AcceptorActivityPlotsWithTree_acceptorsFullCleanedUp.tsv")
gtpred.reactions$source = "GT-Predict"
gtpred.reactions$rate = NA

gtpred.ext = fread("GT-Predict/extensions.AA.tsv", select=c("Acceptor", "Protein", "Reactive"),
                   col.names=c("acceptor", "enzyme", "reaction"))
gtpred.ext$source = "GT-Predict extensions"
gtpred.ext$rate = NA


# from lit we assume all reported rates means there is reactivity since there 
# is too few observations to fit a distribution.
lit = fread("lit/lit.tsv", select=c("Acceptor", "Protein", "Reactive", "KcatPerSec", "Reference"),
            col.names=c("acceptor", "enzyme", "reaction", "rate", "source"))
lit[!is.na(rate), reaction:=1]

hts_ugt.rates = fread("Fatemeh_eval/all-experimental.tsv")
setnames(hts_ugt.rates, "V1", "enzyme")
hts_ugt.rates = melt(hts_ugt.rates, "enzyme", variable.name="acceptor", value.name="rate")
hts_ugt.rates$source = "HTS_UGT"
tmh.rates = fread("inhouse/pTMH.tsv", drop=c("sequence", "species"))
tmh.rates$source = "pTMH"
rates = rbind(hts_ugt.rates, tmh.rates)
rates$reaction = NA
# based on work in rate2bool results folder and discussion with David we 
# decided to set the threshold for reactivity per enzyme 
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

gtpred.reactions$NormEnz = NA
gtpred.ext$NormEnz = NA
lit$NormEnz = NA
DT = rbind(rates, gtpred.reactions, gtpred.ext, lit)

# add CIDs
nBefore = nrow(DT)
# Currently discarding entries without CID, that does have SMILES
raw2cid = na.omit(fread("reactions/rawAcceptor_cid_title.tsv", sep='\t', drop="smiles"))
# Losing these:
# DT[! acceptor %in% raw2cid$raw, as.character(unique(acceptor))]
DT = merge(DT, raw2cid, by.x="acceptor", by.y="raw")
message("CID annotation of reactions: ", nBefore, " -> ", nrow(DT))
DT[, acceptor:=NULL]
setnames(DT, "title", "acceptor")

# Is there conflict between rate data and GT-Predict data and
# what is the lowest rate for an enzyme where a postive datapoint is found from GT-Predict?
rateVreact = merge(DT[!is.na(rate) &  is.na(reaction), .SD, .SDcols=!c("reaction", "acceptor")],
                   DT[ is.na(rate) & !is.na(reaction), .SD, .SDcols=!c("rate", "NormEnz")],
                   by=c("cid", "enzyme"), suffixes=c(".rate", ".reaction"))
rateVreact[order(enzyme, rate)]

library(ggplot2)

ggplot(rateVreact, aes(x=rate, y=reaction, shape=acceptor, color=as.factor(NormEnz))) +
    facet_wrap(enzyme ~ .) +
    geom_point(size=5) +
    scale_shape_manual(values=c(0, 1, 2, 3, 4, 5, 6, 20, 18), name="Acceptor") +
    scale_y_discrete(name="GT-Predict data", labels=c("unreactive", "reactive")) +
    scale_x_continuous(name="Dataset 1 rate") +
    scale_color_discrete(name="Classification based on rate", labels=c("unreactive", "unclear", "reactive"))
ggsave("reactions/GTpredVrate.pdf", height=4, width=9)

# looks actually pretty consistent. We keep data from both sides in spite of 5 conflicting points, 
# since using both will have an effect of somewhere in the middle for a randomforest 
# with its random sampling etc.
DT[!is.na(NormEnz) & is.na(reaction), reaction:=NormEnz]
DT[,NormEnz:=NULL]

sebastian = fread("lit/reactions.tsv")
sebastian$source = "sebastian"
# we define lit data where a rate is reported as reactive per discussion with David
sebastian[!is.na(rate), reaction:=1]
# make sure there are no spaces written around CIDs
sebastian[, cid:=as.integer(cid)]

DT = rbind(DT, sebastian)
DT[reaction=="yes", reaction:=1]
DT[reaction=="no", reaction:=0]
DT = unique(DT)

# NOTE that the CIDs are not validated yet, i.e. may be for chemicals that contain CO-. See results/*validat*
fwrite(DT, "reactions/reactions.tsv", sep='\t')

