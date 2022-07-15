#!/usr/bin/env Rscript
library(data.table)
library(ggplot2)
library(fdrtool)
suppressPackageStartupMessages(library(here))

dt = fread(paste0(here(), "/data/reactions/reactions.tsv"))
dt = dt[!is.na(rate) & source%in%c("pTMH", "HTS_UGT")]


fdrp = function(x, statistic) fdrtool(x, statistic, plot=F, verbose=F)$pval
fdrq = function(x, statistic="pvalue") fdrtool(x, statistic, plot=F, verbose=F)$qval


# vals = dt[enzyme=="At_71C1", rate]
vals = dt[, rate]
vals = sort(vals, decreasing=T)

plot(vals, type='l')

signif = vals[fdrtool(vals, plot=F)$pval < 0.05]
points(1:length(signif), signif)
signif = vals[fdrtool(vals, plot=F)$qval < 0.05]
points(1:length(signif), signif)
signif = vals[fdrtool(vals, plot=F)$lfdr < 1]

cutoff = fdrtool(vals, plot=F)$param[1]
lines(c(0, 2000), c(cutoff, cutoff))


signif = vals[pnorm(vals, mean=0, sd=mean(vals), lower.tail=F) < 0.05]
points(1:length(signif), signif)
signif = vals[phalfnorm(vals, theta=sd2theta(mean(vals)), lower.tail=F) < 0.05]
points(1:length(signif), signif)



dt = dt[, .(enzyme, acceptor, rate)]


dt[,pHalfnormEnz:=phalfnorm(rate, theta=sd2theta(mean(rate)), lower.tail=F), by=enzyme]
dt[,qHalfnormEnz:=p.adjust(pHalfnormEnz), by=enzyme]
dt[,pNormEnz:=pnorm(rate, sd=mean(rate), lower.tail=F), by=enzyme]
dt[,qNormEnz:=p.adjust(pNormEnz), by=enzyme]
# dt[,pFdrtoolEnz:=fdrp(rate, "normal"), by=enzyme]
# dt[,qFdrtoolEnz:=fdrq(rate, "normal"), by=enzyme]
# dt[is.na(dt)] = 1
dt[(pHalfnormEnz < 0.05) & (qHalfnormEnz < 0.05), HalfnormEnz:=T]
dt[(pHalfnormEnz > 0.05) & (qHalfnormEnz > 0.05), HalfnormEnz:=F]
dt[(pNormEnz < 0.05) & (qNormEnz < 0.05), NormEnz:=T]
dt[(pNormEnz > 0.05) & (qNormEnz > 0.05), NormEnz:=F]
# dt[(pFdrtoolEnz < 0.05) & (qFdrtoolEnz < 0.05), FdrtoolEnz:=T]
# dt[(pFdrtoolEnz > 0.05) & (qFdrtoolEnz > 0.05), FdrtoolEnz:=F]
dt[is.na(dt)] = 0.5

# try stringent vs non stringent (p and q) and discard unsure datapoints between the two thresholds

dt.melt = melt(dt, c("enzyme", "acceptor", "rate"), variable.name="method", value.name="signif")
dt.melt[startsWith(as.character(method), "p"), signif:=signif < 0.05]
dt.melt[startsWith(as.character(method), "q"), signif:=signif < 0.05]
dt.melt[, signif:=as.factor(signif)]

# prettier plotting labels
dt.melt[signif==0,   signif:="neither"]
dt.melt[signif==1,   signif:="both p and q"]
dt.melt[signif==0.5, signif:="only p"]

plot_rank = function(meth) {
    ggplot(dt.melt[method==meth], aes(rank(rate), rate, color=signif, shape=signif)) +
        facet_wrap(~enzyme) +
        geom_point(alpha=0.5) +
        theme_linedraw() +
        ggtitle(meth) +
        theme(panel.grid.major=element_line(color="gray"),
            panel.grid.minor=element_line(color="lightgray")) +
        scale_color_manual(name="Significant", values=c("grey", "red", "blue")) +
        scale_shape_manual(name="Significant", values=c(1, 4, 3))
}

plot_rank("HalfnormEnz")
ggsave("thres_HalfnormEnz.pdf", width=9, height=9)
plot_rank("NormEnz")
ggsave("thres_NormEnz.pdf", width=9, height=9)
# plot_rank("FdrtoolEnz")
# ggsave("thres_FdrtoolEnz.pdf", width=9, height=9)


dt.melt[(method=="HalfnormEnz") & (signif==1), .N]
dt.melt[(method=="HalfnormEnz") & (signif!=0), .N]
dt.melt[(method=="NormEnz") & (signif==1), .N]
dt.melt[(method=="NormEnz") & (signif!=0), .N]
dt.melt[(method=="NormEnz") & (signif==0), .N]
dt.melt[(method=="NormEnz") & (signif==0.5), .N]
dt.melt[(method=="FdrtoolEnz") & (signif==1), .N]
dt.melt[(method=="FdrtoolEnz") & (signif!=0), .N]


ggplot(dt.melt[method=="NormEnz"],
    aes(rank(rate), rate, color=signif, shape=signif)) +
    facet_wrap(~enzyme, scales="free") +
    geom_point() +
    theme_linedraw() +
    theme(panel.grid.major=element_line(color="gray"),
          panel.grid.minor=element_line(color="lightgray")) +
    scale_color_manual(name="Significant", values=c("grey", "red", "blue")) +
    scale_shape_manual(name="Significant", values=c(1, 4, 3))

ggsave("rate2bool_all.pdf", width=9, height=9)

ggplot(dt.melt[method=="NormEnz" & enzyme%in%c("Dc_71F5", "Mt_78G1", "Rc_GT1")],
    aes(rank(rate), rate, color=signif, shape=signif)) +
    facet_wrap(~enzyme, scales="free") +
    geom_point() +
    theme_linedraw() +
    theme(panel.grid.major=element_line(color="gray"),
          panel.grid.minor=element_line(color="lightgray")) +
    scale_color_manual(name="Significant", values=c("grey", "red", "blue")) +
    scale_shape_manual(name="Significant", values=c(1, 4, 3))

ggsave("rate2bool_3.pdf", width=9, height=2.8)


