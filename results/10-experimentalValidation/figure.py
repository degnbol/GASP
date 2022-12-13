#!/usr/bin/env python3
import pandas as pd
from plotly import express as px

df = pd.read_table("pred_seqs_blosum62Amb-yield.tsv")

px.scatter(df, x="Yield", y="pred", symbol="Description", color="Polyphenol Name")


