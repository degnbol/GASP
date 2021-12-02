#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(webchem))
suppressPackageStartupMessages(library(here))
setwd(paste0(here(), "/data"))

# read
hts_ugt.raw = colnames(fread("../Fatemeh_eval/all-experimental.tsv"))[-1]
tmh.raw = unique(fread("pTMH.tsv", select="acceptor")$acceptor)
gtpred.raw = unique(fread("gtpred_reactions.tsv", select="acceptor")$acceptor)


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


# manually found after doing the get_cid. Put up here before others to save time, since get_cid is a slow call.
raw2cid = rbind(
    c("Gingerol", 442793),
    c("AbscisicAcid", 5280896),
    c("Abscisic acid", 5280896),
    c("Pinoresinol", 73399),
    c("Borneol", 64685),
    c("Eucamalol", 12426239),
    c("Ginkgolide B", 65243),
    c("Hinokiol", 12310492),
    c("Homoharringtonine", 285033),
    c("Huperzine A", 854026),
    c("Ilicic acid", 11876195),
    c("Ingenol-3-Angelate", 23581946),
    c("Isosteviol", 99514),
    c("Madecassic acid", 73412),
    c("Naringenin", 932),
    c("Nerolidol", 5284507),
    c("Pseudolaric acid B", 53377390),
    c("Theaflavin", 135403798),
    c("DihydroZeatin", 32021),
    c("GDPFuc", 135402013),
    c("GDPMan", 135398627),
    c("Gibberellin A4", 92109),
    c("Indole 3-acetate", 801),
    c("Gibberellin A3", 6466),
    c("UDPGlc", 8629),
    c("2-Methyl-umbelliferone", 5280567),  # 2-methyl doesn't seem to exist and it says 4-methyl in the GT Predict paper so it is a typo
    c("()- cis, trans Abscisic acid", 5280896),  # looked in the GT Predict paper
    c("a-cyano-4-hydroxyl-cinnamic acid", 5328791),
    c("Trans-Zentin-Glucose", 449093), # from GT-Predict. They only refer to Trans-Zeatin in paper.
    c("MUGlcNAc", 2733787), # or 118328, but they have the same canonical SMILES
    c("BocCysThrOMe", NA),
    c("1-Thio-S-cyanomethyl-N-acetyl-D-glucosamine", NA),
    c("?-GlcOBn", NA),
    c("a-ManOBn", NA),
    c("a-ManOCH2Bn", NA),
    c("a-ManOPh", NA),
    c("a-ManOPMP", NA),
    c("a-ManOBn(pNO2)", NA),
    c("a-ManOPhF5", NA),
    c("a-ManOBnF5", NA),
    c("ManSTol", NA), # mannose-S-toluene ?
    c("UDP5SGlc", 8629),  # I think based on looking at the GT-Predict paper that this is just UDP-glucose bound at a 5S location 
    c("dTDPXyl", 122707150), # best match but doesn't seem like a perfect match?
    c("UDPGlcNAc", 445675),
    c("dTDPGlc", 443210),
    c("GDPGlc", 135398625),
    c("UDPMan", 448873),
    c("UDPRha", 49852436),
    c("dTDPRha", 49852346),
    c("FuranThiol", 143754), # best match but not perfect
    c("Lovastatin,Terpineol", NA) # which or both
)
colnames(raw2cid) = c("raw", "cid")


# try to standardize names

accs.raw = unique(c(hts_ugt.raw, tmh.raw, gtpred.raw))
acceptors = data.table(raw=accs.raw, processed=NA)

acceptors2 = data.table(rbind(
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
colnames(acceptors2) = c("raw", "processed")

acceptors = rbind(
    acceptors,
    acceptors2,
    acceptors[grep(" (", fixed=T, raw), .(raw, processed=sub(".* \\((.*)\\)", "\\1", raw))],
    acceptors[grep(" (", fixed=T, raw), .(raw, processed=sub(" \\(.*", "", raw))],
    acceptors[, .(raw, processed=subs(c("^ ", "?-", "(?)-", "()-"), "", fixed=c(F,T,T,T), raw))],
    acceptors[grep("[a-z][A-Z]$", raw), .(raw, processed=sub("A$", " A", raw))],
    acceptors[grep("^GDP[A-Z]", raw), .(raw, processed=sub("^GDP", "GDP-", raw))]
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
    c("a-", "alpha-")
)

acceptors = rbind(acceptors, acceptors[, .(raw, processed=subs(substitutions[,1], substitutions[,2], processed))])
acceptors[, processed:=sub("--", "-", processed)]
acceptors[, processed:=sub("-$", "", processed)]
acceptors = acceptors[!processed%in%c("DCA")]
acceptors = unique(na.omit(acceptors))

cids = get_cid(acceptors[!raw%in%raw2cid[,"raw"], processed], match="na")
setnames(cids, "query", "processed")

acceptors = unique(merge(acceptors, na.omit(cids))[, cid, by=raw])
# remove any raw2cid mapping that has ambiguity
acceptors = acceptors[!raw%in%acceptors[,.N, by=raw][N>1, raw]]

# add anything else manually to be sure it's OK, i.e. anything not found or found with multiple CIDs. 
accs.raw[!accs.raw%in%acceptors$raw]

acceptors = na.omit(unique(rbind(acceptors, raw2cid)))

fwrite(acceptors, "raw-acceptor2cid.tsv", sep='\t')

