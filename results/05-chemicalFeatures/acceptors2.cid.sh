#!/usr/bin/env zsh
# The first acceptors were training data curated from public sources and from two inhouce datasets.
# Then experiments were carried out to assess the performance with the data named 20220215.
# CIDs from both are combined here for any testing newer than these.
cat acceptors.cid 20220215.cid | sort -nu > acceptors2.cid
# manually add acceptors without pubchem ID
echo 'O=C(N[C@@H](CS)C(N[C@H](C(OC)=O)[C@@H](O)C)=O)OC(C)(C)C
OC[C@H]1O[C@H](OC2=C(F)C(F)=C(F)C(F)=C2F)[C@@H](O)[C@@H](O)[C@@H]1O
OC[C@H]1O[C@H](OCC2=C(F)C(F)=C(F)C(F)=C2F)[C@@H](O)[C@@H](O)[C@@H]1O
OC[C@H]1O[C@@H](SCCCCCCCCCCCCCCC(O)=O)[C@H](NC(C)=O)[C@@H](O)[C@@H]1O
OC[C@H]1O[C@@H](O[C@H]2[C@H](O)[C@@H](NC(C)=O)[C@H](SCC#N)O[C@@H]2CO)[C@H](O)[C@@H](O)[C@H]1O
O[C@H]1[C@H](O[C@H]2[C@@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O2)[C@H](O)[C@@H](OC3=CC=CC=C3)O[C@@H]1CO[C@@H]4[C@@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O4' \
    >> acceptors2.cid
