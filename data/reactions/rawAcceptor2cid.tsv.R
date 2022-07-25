#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(webchem))
suppressPackageStartupMessages(library(here))
setwd(paste0(here(), "/data"))

# functions
subs = function(patterns, replacements, x, fixeds=F) {
    n = length(patterns)
    if(length(fixeds) == 1) fixeds = rep(fixeds, n)
    if(length(replacements) == 1) replacements = rep(replacements, n)
    o = x
    for(i in 1:n) {
        o = sub(patterns[i], replacements[i], o, fixed=fixeds[i])
    }
    o
}

# read files with raw names that needs converting to cid
hts_ugt.raw = colnames(fread("Fatemeh_eval/all-experimental.tsv"))[-1]
tmh.raw = unique(fread("inhouse/pTMH.tsv", select="acceptor")$acceptor)
gtpred.raw = unique(fread("reactions/gtpred_reactions.tsv", select="acceptor")$acceptor)
gtpred_ext.raw = unique(fread("GT-Predict/extensions.tsv"))
lit.raw = unique(fread("lit/lit.tsv")$Acceptor)

accs.raw = unique(c(hts_ugt.raw, tmh.raw, gtpred.raw, lit.raw))

# manual effort
raw2cid.manual = fread("reactions/rawAcceptor2cid_manual.tsv", drop=c("smiles", "comment"))

# save curation time by only running webchem for acceptors that haven't already been looked up.
# Do this by reading the result from this script if present and skipping those entries in the lookup.
if(file.exists("reactions/rawAcceptor2cid.tsv")) {
    raw2cid = fread("reactions/rawAcceptor2cid.tsv")
    raw2cid = rbind(raw2cid, raw2cid.manual)
} else {
    raw2cid = raw2cid.manual
}

# define alternative spellings etc for queries corresponding to each raw acceptor string.

raw2query = data.table(raw=accs.raw, query=NA)

raw2query.manual = data.table(rbind(
    c("BenzAdenine", "Benzyladenine"),
    c("DiMeBenzAcid", "Dimethylbenzoic acid"),
    c("CinAcid", "Cinnamic acid"),
    c("CinAlc", "Cinnamyl alcohol"),
    c("ConAlc", "Coniferyl alcohol"),
    c("DCA", "3,4-dichloroaniline"),
    c("DCP", "3,4-dichlorophenol"),
    c("DCT", "3,4 dichlorothiophenol"),
    c("TCP", "Trichlorophenol"),
    c("Gibberellin A3", "Gibberellic acid"),
    c("HomovanilAcid", "Homovanillic acid"),
    c("3,4-dihyroxy cinnamic acid", "3,4-dihydroxycinnamic acid"),
    c("4-hydroxy 3-methoxy cinnamic acid", "4-hydroxy-3-methoxycinnamic acid"),
    c("Coumaryl alcohol", "p-Coumaryl alcohol")
))
colnames(raw2query.manual) = c("raw", "query")

# for each raw acceptor string we get multiple acceptor queries that can be fed to pubchem,
# e.g. "a (b)" means both a and b are search strings.
raw2query = rbind(
    raw2query,
    raw2query.manual,
    raw2query[grep(" (", fixed=T, raw), .(raw, query=sub(" \\(.*", "", raw))],
    raw2query[, .(raw, query=subs(c("^ ", "?-", "(?)-", "()-"), "", fixed=c(F,T,T,T), raw))],
    raw2query[grep("[a-z][A-Z]$", raw), .(raw, query=sub("A$", " A", raw))],
    raw2query[grep("^GDP[A-Z]", raw), .(raw, query=sub("^GDP", "GDP-", raw))]
)

substitutions = rbind(
    c("Acid$", " acid"),
    c("Alc$", " alcohol"),
    c("^Glc", "Glucose-"),
    c("Glc", "-Glucose-"),
    c("Xyl$", "-xylose"),
    c("Fuc$", "-fucose"),
    c("Me$", "Methyl"),
    c("Man([A-Z])", "mannose\\1"),
    c("Me([A-Z])", "Methyl\\1"),
    c("hydroxylbenzoic", "hydroxybenzoic"),
    c("hydroxyl-", "hydroxy"),
    c("hydroxy ", "hydroxy"),
    c("coumerin", "coumarin"),  # typo
    c("Zentin", "zeatin"),   # typo
    c(" ([0-9])", "-\\1"),
    c("a-", "alpha-"),
    c("α", "alpha"),
    c("β", "beta"),
    c("hydroximate$", "hydroximic acid"), # conjugate acid
    c("benzoate$", "benzoic acid") # conjugate acid
)

raw2query = rbind(raw2query, raw2query[, .(raw, query=subs(substitutions[,1], substitutions[,2], query))])
raw2query[, query:=sub("--", "-", query)]
raw2query[, query:=sub("-$", "", query)]
raw2query = raw2query[!query%in%c("DCA")]
raw2query = unique(na.omit(raw2query))

# skip queries for raw acceptor names that have already been curated.
queries = raw2query[! raw %in% raw2cid$raw]

nRawLeft = length(unique(queries$raw))
message("Running ", nrow(queries), " get_cid queries from ", nRawLeft, " raw chemical names...")
# match="na" means we don't get a match if there are multiple returned by 
# pubchemwe don't get a match if there are multiple returned by pubchem.
cids = get_cid(queries$query, match="na")
cids = na.omit(cids)
raw.matched = unique(merge(cids, queries, all.x=T))$raw
message(length(raw.matched), "/", nRawLeft, " matched raw acceptor names.")
queries = queries[! raw %in% raw.matched]

# second round of matching is more flexible, in cases where there was no match for the more strict search.
# Here, add the contents of parenthesis in raw string:
raw2query.paren = raw2query[grep(" (", fixed=T, raw), .(raw, query=sub(".* \\((.*)\\)", "\\1", raw))]
raw2query = rbind(raw2query, raw2query.paren)
queries = rbind(queries, raw2query.paren)
                  
# in cases where there are multiple matches we can take the first match,
# but only if we didn't manage to find a match with any of the variants 
# ("query") that describes a given raw string.
nRawLeft = length(unique(queries$raw))
message("Running ", nrow(queries), " get_cid flexible queries from ", nRawLeft, " raw chemical names...")
cids.flex = get_cid(queries$query, match="first")
cids.flex = na.omit(cids.flex)
raw.matched = unique(merge(cids.flex, queries, all.x=T))$raw
message(length(raw.matched), "/", nRawLeft, " matched raw acceptor names.")
queries = queries[! raw %in% raw.matched]

# remaining names have to be added manually
lacking = unique(queries$raw)
message(length(lacking), " raw strings not matched:")
message(paste(lacking, collapse='\n'))

cids = merge(cids, raw2query, by="query")
raw2cid.new = unique(data.table(merge(cids, raw2query))[, cid, by=raw])
# multiple matches:
multiMatch = raw2cid.new[raw%in%raw2cid.new[,.N, by=raw][N>1, raw]]
nMultiMatch = length(unique(multiMatch$raw))
if(nMultiMatch > 0) message(nMultiMatch, " raw strings with multiple matches.")
# remove them
raw2cid.new = raw2cid.new[!raw%in%raw2cid.new[,.N, by=raw][N>1, raw]]

raw2cid = na.omit(unique(rbind(raw2cid, raw2cid.new)))
fwrite(raw2cid, "reactions/rawAcceptor2cid.tsv", sep='\t')

