#!/bin/bash
#SBATCH --job-name=featSel
#SBATCH --partition=mig
#SBATCH -c 128
#SBATCH --mem=64G
#SBATCH --time=12-00:00:00
#SBATCH -o feature_selection.out
#SBATCH -e feature_selection.err
#SBATCH --mail-user=cdmadsen@student.unimelb.edu.au
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80

if [ -z "$SLURM_JOB_ID" ]; then
    me=`basename $0`
    echo "Submit with sbatch $me"
    exit 1
fi

julia -t 128 ./feature_selection.tsv.gz.jl

my-job-stats -a -n -s
