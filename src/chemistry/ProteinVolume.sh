#!/usr/bin/env zsh
# USAGE: src/chemistry/ProteinVolume.sh "pdbs" > $WORK/volumes.tsv
# volumes.tsv will have columns Protein,SolventExcludedVolume,VanDerWaalsVolume
# use PDBs in folder given in only argument to make volumes and also writes intermediary logs to pdbs/OutputDir_pdbs_*.txt
# we can get solvent-excluded volume with ProteinVolume https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0531-2
# http://gmlab.bio.rpi.edu/download.php
WORK=$1

ROOT=`git root`
TOOL="$ROOT/tools/ProteinVolume_1.3" # dir for ProteinVolume

# rm previous work
# (N) to glob without error if no file is found
pdb_logs=($1/OutputDir_pdbs_*.txt(N))
[ "$pdb_logs" ] && rm $pdb_logs

# number of structures
N=`ls $WORK/*.pdb | wc -l | xargs`

java -jar $TOOL/ProteinVolume_1.3.jar -het --radiusFileName $TOOL/bondi.rad $WORK |
    tee ProteinVolume.tmp.tsv | $ROOT/src/progress.sh "$N" 6 1>&2 # STDERR

sed '1,6d' ProteinVolume.tmp.tsv | sed $'s/   */\t/g' |
    mlr --tsv rename -g -r ' ,' then cut -x -f 'TimeTaken(ms)' then \
    rename 'TotalVolume(A3),SolventExcludedVolume,VDWVolume,VanDerWaalsVolume,Protein,cid' then \
    sort -n cid | tr ',' '.' # make sure it uses . as decimal then STDOUT

# cleanup
rm ProteinVolume.tmp.tsv
