#!/usr/bin/env zsh
# find the bases that needs correcting.
./CO-.tsv.sh
echo "Manually look up the conjugates and write conjugates.tsv"
# correct reactions cids from base->acid using conjugates
./reactions.tsv.R
# get corrected acceptors
./acceptors.tsv.sh
# assert that all is corrected
if [ `./getCO-.py < acceptors.tsv | sed 1d` ]; then
    echo "Acceptors not converted from conjugate base to acid:"
    ./getCO-.py < acceptors.tsv
    exit 1
fi
    

