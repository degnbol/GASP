#!/bin/bash
#SBATCH --job-name=featSel
#SBATCH --partition=mig
#SBATCH -c 24
#SBATCH --mem=100G
#SBATCH --time=10:00:00
#SBATCH -o feature_selection.out
#SBATCH -e feature_selection.err
#SBATCH --mail-user=cdmadsen@student.unimelb.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT_80

if [ -z "$SLURM_JOB_ID" ]; then
    me=`basename $0`
    echo "Submit with sbatch $me"
    exit 1
fi

julia -t 24 ./feature_selection.tsv.gz.jl

my-job-stats -a -n -s
