#!/usr/bin/env Rscript
library(data.table)
library(ggplot2)
suppressPackageStartupMessages(library(here))

dt = fread(paste0(here(), "data/reactions/reactions.tsv"))
dt = dt[!is.na(rate)]
dt = dt[cid%in%c(7257, 27582, 31553, 5280343)]

unique(dt[, .N, by=enzyme]$N)
dt[acceptor=="quercetin", acceptor:="Quercetin"]

ggplot(dt, aes(x=enzyme, y=rate, fill=acceptor)) +
    geom_col(position="dodge") +
    theme_linedraw() +
    theme(axis.text.x=element_text(angle=90, vjust=0.5),
          panel.grid.major=element_line(color="gray"), panel.grid.minor=element_line(color="lightgray"))

ggsave("rates.pdf", width=10, height=5)


