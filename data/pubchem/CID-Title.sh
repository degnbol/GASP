#!/usr/bin/env zsh
path=`git root`/data/pubchem/CID-Title.gz
wget -O $path ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-Title.gz
gunzip -k $path
