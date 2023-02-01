#!/bin/bash
#SBATCH --job-name=featSel2
#SBATCH --partition=mig
#SBATCH -c 72
#SBATCH --mem=400G
#SBATCH --time=30-00:00:00
#SBATCH -o feature_selection-2.out
#SBATCH -e feature_selection-2.err
#SBATCH --mail-user=cdmadsen@student.unimelb.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT_80

if [ -z "$SLURM_JOB_ID" ]; then
    me=`basename $0`
    echo "Submit with sbatch $me"
    exit 1
fi

julia -t 72 ./feature_selection-2.tsv.gz.jl

my-job-stats -a -n -s
