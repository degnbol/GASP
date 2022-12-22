#!/usr/bin/env Rscript
if (!require("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("ggtree", force=TRUE)

library("data.table")
library("ggplot2")
suppressPackageStartupMessages(library("ggtree"))
library(patchwork) # simply combine figs with +

abbrev = function(x) {
    stringr::str_trunc(x, 20)
}


tree = read.tree("dataset1.tsv.tree")

plt = ggtree(tree, branch.length="none") + geom_tiplab()
plt = ggtree(tree) + geom_tiplab(align=TRUE, size=0) +
    theme(plot.margin=unit(c(-.5,-1.5,-.5,-1.5), "lines"))

# extract the eznyme ordering from the tree
enzyme_hclust = data.table(plt[["data"]])[isTip==TRUE,][order(y), label]

dt = fread("dataset1.tsv")
acceptor_hclust = readLines("acceptor_hclust.txt")
dt[, acceptor:=factor(acceptor, levels=acceptor_hclust)]
dt[, enzyme:=factor(enzyme, levels=enzyme_hclust)]

# have NA entries for enzyme-acceptor pairs not recorded to control its coloring
allVall = expand.grid(acceptor=unique(dt$acceptor), enzyme=unique(dt$enzyme))
dt = merge(dt, allVall, all=TRUE)

dt[reaction==0.5, reactive:="inconclusive"]
dt[reaction==0, reactive:="no"]
dt[reaction==1, reactive:="yes"]
dt[is.na(reaction), reactive:="unmeasured"]
dt[, reactive:=factor(reactive, levels=c("unmeasured", "no", "inconclusive", "yes"))]

dt[!is.na(rate), label:=sprintf("%.2f", rate)]
dt[rate <= 0.005, label:="0"]

dt[,acceptor_trunc:=stringr::str_trunc(as.character(acceptor), 30)]

dt[, color:=rate]
dt[rate<= 0.005, color:=0]
# dt[is.na(color), color:=0]
# avoid na in log(color)
dt[color==0, color:=min(dt[color>0,color])]

# a few options for colors
fill_vals = c("white", "red", "blue", "green");    color_low = "lightgray"; color_high = "black"
fill_vals = c("white", "blue", "black", "yellow"); color_low = "lightgray"; color_high = "black"
fill_vals = c("white", "white", "gray", "black");  color_low = "red";       color_high = "green"

plt2 = ggplot(dt, aes(x=acceptor, y=enzyme, fill=reactive, color=color, label=label)) +
    geom_tile(linewidth=0.5) +
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
          # no change from vjust
          # axis.text.x=element_text(angle=90, hjust=0, vjust="inward"),
          axis.text.x=element_blank(),
          axis.text.y=element_text(hjust=0, margin=margin(c(-1,-1,-1,-1)), size=7),
          axis.ticks=element_blank(), axis.title=element_blank(),
          plot.margin = unit(c(7.1, 0, 0, 0), "lines"), # make space for drawing acceptor names outside plotting area
          legend.key.size=unit(1,"lines"), legend.title=element_text(size=10), legend.text=element_text(size=8)) +
    scale_x_discrete(expand=c(0,0), position="top", label=abbrev) +
    scale_y_discrete(expand=c(0,0)) +
    scale_fill_manual(values=fill_vals) +
    scale_color_gradient(name="rate", labels=function(x) sprintf("%.3f", x), low=color_low, high=color_high, trans="log", na.value="lightgray") +
    geom_text(size=2, angle=90, fontface="bold") +
    # hack solution since vjust doesn't work
    geom_text(mapping=aes(x=acceptor, label=acceptor_trunc), inherit.aes=FALSE, size=2.5, y=length(levels(dt$enzyme))+0.5+0.1, angle=90, hjust=0) +
    coord_cartesian(clip='off') # draw outside plot


plt + plt2 + plot_layout(widths=c(1,4))

ggsave("dataset1_3.pdf", width=11, height=6.2)















