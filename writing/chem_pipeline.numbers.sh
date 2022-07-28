#!/usr/bin/env zsh
# print some numbers to put in figure.
infile=`git root`/data/reactions/reactions.tsv

gtpredFilt='$source == "GT-Predict" || $source == "GT-Predict extensions"'
litFilt='$source != "GT-Predict" && $source != "GT-Predict extensions" && $source != "pTMH" && $source != "HTS_UGT"'

echo "GT-Predict"
echo -n "  reactions = "
mlr -t --from $infile filter $gtpredFilt + count | sed 1d
echo -n "  acceptors = "
mlr -t --from $infile filter $gtpredFilt + uniq -f cid + count | sed 1d
echo -n "  enzymes   = "
mlr -t --from $infile filter $gtpredFilt + uniq -f enzyme + count | sed 1d

echo "HTS_UGT"
echo -n "  reactions = "
mlr -t --from $infile filter '$source == "HTS_UGT"' + count | sed 1d
echo -n "  acceptors = "
mlr -t --from $infile filter '$source == "HTS_UGT"' + uniq -f cid + count | sed 1d
echo -n "  enzymes   = "
mlr -t --from $infile filter '$source == "HTS_UGT"' + uniq -f enzyme + count | sed 1d

echo "pTMH"
echo -n "  reactions = "
mlr -t --from $infile filter '$source == "pTMH"' + count | sed 1d
echo -n "  acceptors = "
mlr -t --from $infile filter '$source == "pTMH"' + uniq -f cid + count | sed 1d
echo -n "  enzymes   = "
mlr -t --from $infile filter '$source == "pTMH"' + uniq -f enzyme + count | sed 1d

echo "lit"
echo -n "  reactions = "
mlr -t --from $infile filter $litFilt + count | sed 1d
echo -n "  acceptors = "
mlr -t --from $infile filter $litFilt + uniq -f cid + count | sed 1d
echo -n "  enzymes   = "
mlr -t --from $infile filter $litFilt + uniq -f enzyme + count | sed 1d

