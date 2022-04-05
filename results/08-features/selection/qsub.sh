qsub -d $PWD -l nodes=1:ppn=30,walltime=80:00:00,mem=8g $PWD/feature_importance.sh
