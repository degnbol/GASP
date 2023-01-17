#!/bin/bash
#SBATCH --job-name=featSel
#SBATCH --partition=physical
#SBATCH -c 128
#SBATCH --mem=128G
#SBATCH --time=20-00:00:00
#SBATCH -o feature_selection.out
#SBATCH -e feature_selection.err

if [ -z "$SLURM_JOB_ID" ]; then
    me=`basename $0`
    echo "Submit with sbatch $me"
    exit 1
fi

julia -t 128 ./feature_selection.tsv.gz.jl

my-job-stats -a -n -s
