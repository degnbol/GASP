#!/usr/bin/env zsh
# All-vs-all Smith-Waterman alignments.
# REQUIREMENTS: emboss water. Install with brew install emboss.
# USAGE: water.sh TARGETS.fa [options] < QUERIES.fa > OUTFILE.tsv
cat - | while read header; read seq
do
    echo "$header\n$seq" |
        water -auto -filter -bsequence ${=@}
done |
    grep -E '^# 1|^# 2|^# Length|^# Identity|^# Gaps|^# Score|^$' |
    sed 's/\/.*//' | tr -d '# (%)' |
    mlr --x2t --ips colon rename '1,Query,2,Target,Identity,Matches' +\
    put '$Identity = $Matches/$Length'
