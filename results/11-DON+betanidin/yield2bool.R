#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
library(ggplot2)
library(fdrtool)
suppressPackageStartupMessages(library(here))

cwd = function(s) paste(here(), s, sep="/")

dt = fread(cwd("results/10-experimentalValidation/experimentYield.tsv"))
yield = dt$Yield

mean(yield)
median(yield)

plot(rank(yield, ties.method="first"), yield, pch='.')


dt[,pNorm:=phalfnorm(Yield, theta=sd2theta(mean(Yield)), lower.tail=F)]
dt[,qNorm:=p.adjust(pNorm)]
dt[(pNorm < 0.05) & (qNorm < 0.05), Significant:="both p and q"]
dt[(pNorm < 0.05) & (qNorm > 0.05), Significant:="only p"]
dt[(pNorm > 0.05) & (qNorm > 0.05), Significant:="neither"]

ggplot(dt, aes(x=rank(Yield, ties.method="first"), y=Yield, color=Significant)) +
    geom_point()

# the fit norm p and q method doesn't look good. trying thres >25 and <12.5 as yes and no in train.tsv.sh


