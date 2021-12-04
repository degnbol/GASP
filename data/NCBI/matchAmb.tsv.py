#!/usr/bin/env zsh
# add lookup for special amino acid codes B, J, Z, and X which mean AsparticAcid/asparagine, leucine/isoleucine, glutamicAcid/glutamine, and any, respectively.
# https://www.ddbj.nig.ac.jp/ddbj/code-e.html
# this is different from default NCBI encoding, which makes them separate codes (see raw/MATCH) which is a bad approach, at least for any we want to use these for.
from os.path import expanduser
import pandas as pd

dir = expanduser("~/biosustain/data/ncbi/")

match = pd.read_table(dir + "match.tsv")
match["B"] = match["J"] = match["Z"] = -1
match.loc[match["_"].isin(["D", "N"]), "B"] = 1/2
match.loc[match["_"].isin(["I", "L"]), "J"] = 1/2
match.loc[match["_"].isin(["E", "Q"]), "Z"] = 1/2
match["X"] = 1/(len(match) - 1)
match.loc[match["_"] == '*', "X"] = -1

match.to_csv(dir + "matchAmb.tsv", sep='\t', index=False)
