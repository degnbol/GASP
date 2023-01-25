#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(here))

# Check consistency between two potential sources of acceptor interactions from GT-Predict.
# ids are still raw (human readable at the end of this script)

setwd(paste0(here(), "/data/GT-Predict"))

df_full = fread("AcceptorActivityPlotsWithTree_acceptorsFullCleanedUp.tsv")
df_intr = data.table::melt(fread("acceptor_interaction_data.txt.tsv"), "acceptor", variable.name="enzyme", value.name="reaction")
# GT-Predict classification system:
# NA=unreactive
# 0=unreactive
# 1=reactive
# 2=unclear
# 3=missing
df_intr[, reaction:=as.double(reaction)]
df_intr[is.na(reaction), reaction:=0]
df_intr[reaction == 2, reaction:=0.5]
df_intr = df_intr[reaction != 3]

df_cmp = merge(df_full, df_intr, by=c("acceptor", "enzyme"))
cor(df_cmp$reaction.x, df_cmp$reaction.y)
# as expected the files are identical when there is overlap in identification
# we use the largest (df_full) in the reactions.tsv.R


