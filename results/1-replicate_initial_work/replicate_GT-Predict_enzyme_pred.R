#!/usr/bin/env Rscript
library(ggplot2)
library(pROC)
library(data.table)
suppressPackageStartupMessages(library(here))
# go to git root
setwd(here())

water = fread("data/water.tsv")

hts_ugt_rates = fread("data/Fatemeh_eval/all-experimental.tsv")
setnames(hts_ugt_rates, "V1", "enzyme")

gt.pred = fread("data/gtpred_reactions_enz.tsv")
# convert away from their weird classification system where 0/missing = not an enzyme-molecule match, 1= match, 2 or 3 is unknown
gt.pred[is.na(gt.pred)] = F
gt.pred[gt.pred == 1] = T
gt.pred[gt.pred == 2 | gt.pred == 3] = NA

# then between the two projects
hts.ugt_acceptors = names(hts_ugt_rates)[-1]
gt.pred_acceptors = gsub(" ", "", gt.pred$acceptor)
for(acc in hts.ugt_acceptors) {
    # match on the whole phrase
    idx = grep(paste0("^", acc, "$"), gt.pred_acceptors, ignore.case=T)
    if(length(idx) == 1) {
        message("matching acceptor from HTS UGT: ", acc, " to GT-Predict: ", gt.pred$acceptor[[idx]])
        gt.pred[idx, acceptor:=acc]
    } else {
        # match on partial phrase
        idx = grep(acc, gt.pred_acceptors, ignore.case=T)
        if(length(idx) == 1) {
            message("matching acceptor from HTS UGT: ", acc, " to GT-Predict: ", gt.pred$acceptor[[idx]])
            gt.pred[idx, acceptor:=acc]
        }
    }
}

acceptors = hts.ugt_acceptors[hts.ugt_acceptors %in% gt.pred$acceptor]
message("matched on ", paste(acceptors, collapse=", "))

# long formats filtering for matched acceptors
gt.pred = melt(gt.pred[acceptor%in%acceptors], "acceptor", variable.name="enzyme", value.name="class")
hts_ugt_rates = data.table::melt(hts_ugt_rates, "enzyme", variable.name="acceptor", value.name="rate")[acceptor%in%acceptors]

top = water[order(-score), .(first=head(gtpred_enzyme, 1), second=tail(head(gtpred_enzyme, 2), 1)), by=enzyme]
DT = merge(top, gt.pred, by.x="first", by.y="enzyme", all.x=T)
setnames(DT, "class", "first_pred")
DT = merge(DT, gt.pred, by.x=c("second", "acceptor"), by.y=c("enzyme", "acceptor"), all.x=T)
setnames(DT, "class", "second_pred")
DT[, pred:=first_pred]
DT[, gtpred_enzyme:=first]
DT[is.na(first_pred), pred:=second_pred]
DT[is.na(first_pred), gtpred_enzyme:=second]
DT = merge(DT, hts_ugt_rates, by=c("enzyme", "acceptor"))


# plot


ggplot(DT, aes(x=enzyme, y=rate, color=factor(pred), shape=acceptor))+
    geom_point(size=5) +
    theme_bw() + 
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          legend.position='right') +
    theme(axis.text.x = element_text(size=14, angle=90), axis.text.y = element_text(size=14),
          axis.title = element_text(size=18, face="bold")) +
    guides(size = FALSE) + 
    scale_color_manual(values=c('#5ab4ac','#8c510a')) +
    scale_shape_manual(values=seq(0,5))


DT[,cor(rate, pred)]
DT[,auc(roc(pred, rate))]

# then we check if we are overestimating performance because of high identity.

DT = water[DT, on=c("enzyme", "gtpred_enzyme")]


ggplot(DT, aes(x=identity, y=rate, color=factor(pred), shape=acceptor)) +
    geom_point(size=5) +
    theme_bw() + 
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          legend.position='right') +
    theme(axis.text.x = element_text(size=14, angle=90), axis.text.y = element_text(size=14),
          axis.title = element_text(size=18, face="bold")) +
    guides(size = FALSE) + 
    scale_color_manual(values=c('#5ab4ac','#8c510a')) +
    scale_shape_manual(values=seq(0,5))
# we don't see a clear tendency

DT[identity<80, cor(rate, pred)]
DT[identity<80, auc(roc(pred, rate))]

DT[,gmedian:=rate>median(rate)]
DT[,correct:=gmedian==pred]
ggplot(DT, aes(x=identity, fill=correct)) +
    geom_histogram(bins=10)

DT[,cor(correct, identity)]



