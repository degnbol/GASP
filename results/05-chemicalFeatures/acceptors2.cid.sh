#!/usr/bin/env zsh
# The first acceptors were training data curated from public sources and from two inhouce datasets.
# Then experiments were carried out to assess the performance with the data named 20220215.
# CIDs from both are combined here for any testing newer than these.
cat acceptors.cid 20220215.cid | sort -nu > acceptors2.cid
