#!/usr/bin/env zsh
# USAGE: src/chemistry/ProteinVolume.sh "pdbs" > $WORK/volumes.tsv
# volumes.tsv will have columns Protein,SolventExcludedVolume,VanDerWaalsVolume
# use PDBs in folder given in only argument to make volumes and also writes intermediary logs to pdbs/OutputDir_pdbs_*.txt
# we can get solvent-excluded volume with ProteinVolume https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0531-2
# http://gmlab.bio.rpi.edu/download.php
WORK=$1

TOOL="$0:h/../../tools/ProteinVolume_1.3" # dir for ProteinVolume

# rm previous work
# (N) to glob without error if no file is found
pdb_logs=($1/OutputDir_pdbs_*.txt(N))
[ "$pdb_logs" ] && rm $pdb_logs

# number of structures
N=`ls $WORK/*.pdb | wc -l | xargs`

java -jar $TOOL/ProteinVolume_1.3.jar -het --radiusFileName $TOOL/bondi.rad $WORK |
    tee $WORK/ProteinVolume.tmp.tsv | $0:h/../progress.sh "$N" 6 1>&2 # STDERR

# has info at the top of file and may print warnings within data, e.g.
# "Reading hydrogens is turned on, but couldn't find any hydrogens in ..."
# So we grep for lines that have 6 columns separated by spaces, where there the 5 last columns have to be integer or float.
grep -E '^[^ ]+ +([0-9.]+ +){4}' $WORK/ProteinVolume.tmp.tsv | sed $'s/   */\t/g' |
    mlr -t label 'cid,SolventExcludedVolume,VoidVolume,VanDerWaalsVolume,PackingDensity,Time_ms' +\
    sort -n cid | tr ',' '.' # make sure it uses . as decimal then STDOUT

# cleanup
# rm $WORK/ProteinVolume.tmp.tsv
