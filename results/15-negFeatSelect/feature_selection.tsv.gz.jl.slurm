#!/bin/bash
#BATCH --job-name=featSel
#BATCH --partition=mig
#BATCH -c 128
#BATCH --mem=64G
#BATCH --time=12-00:00:00
#BATCH -o feature_selection.out
#BATCH -e feature_selection.err

if [ -z "$SLURM_JOB_ID" ]; then
    me=`basename $0`
    echo "Submit with sbatch $me"
    exit 1
fi

julia -t 128 ./feature_selection.tsv.gz.jl

my-job-stats -a -n -s
