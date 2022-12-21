#!/usr/bin/env Rscript
if (!require("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("ggtree", force=TRUE)

library("data.table")
library("ggplot2")
suppressPackageStartupMessages(library("ggtree"))

abbrev = function(x) {
    stringr::str_trunc(x, 20)
}


tree = read.tree("dataset1.tsv.tree")

plt = ggtree(tree, branch.length="none") + geom_tiplab()
plt = ggtree(tree) + geom_tiplab(align=TRUE)

# extract the eznyme ordering from the tree
enzyme_hclust = data.table(plt[["data"]])[isTip==TRUE,][order(y), label]

dt = fread("dataset1.tsv")
acceptor_hclust = readLines("acceptor_hclust.txt")
dt[, acceptor:=factor(acceptor, levels=acceptor_hclust)]
dt[, enzyme:=factor(enzyme, levels=enzyme_hclust)]
dt[reaction==0.5, reactive:="inconclusive"]
dt[reaction==0, reactive:="no"]
dt[reaction==1, reactive:="yes"]
dt[, reactive:=factor(reactive, levels=c("no", "inconclusive", "yes"))]
dt[!is.na(rate), label:=sprintf("%.2f", rate)]
dt[rate <= 0.005, label:="0"]

dt[, color:=rate]
dt[rate<= 0.005, color:=0]
dt[is.na(color), color:=0]


# have NA entries for enzyme-acceptor pairs not recorded to control its coloring
allVall = expand.grid(acceptor=unique(dt$acceptor), enzyme=unique(dt$enzyme))
dt = merge(dt, allVall, all=TRUE)

ggplot(dt, aes(x=acceptor, y=enzyme, fill=reactive, color=color, label=label)) +
    geom_tile(linewidth=0.5) +
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
          axis.text.x=element_text(angle=90, hjust=0)) +
    scale_x_discrete(expand=c(0,0), position="top", label=abbrev) +
    scale_y_discrete(expand=c(0,0)) +
    scale_fill_discrete(na.value='white') +
    scale_colour_gradient(low="lightgray", high="black", trans="log", na.value="lightgray") +
    geom_text(size=2, angle=90, fontface="bold")

# TODO: combine tree and heatmap, remove y axis title.
