#!/usr/bin/env zsh
# USAGE: src/chemistry/ProteinVolume.sh "pdbs" > volumes.tsv
# volumes.tsv will have columns Protein,SolventExcludedVolume,VanDerWaalsVolume
# use PDBs in folder given in only argument to make volumes and also writes intermediary logs to pdbs/OutputDir_pdbs_*.txt
# we can get solvent-excluded volume with ProteinVolume https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0531-2
# http://gmlab.bio.rpi.edu/download.php
SRC=`git root`"tools/ProteinVolume_1.3" # dir for ProteinVolume

# rm previous work, (N) to glob without error if no file is found
pdb_logs=($1/OutputDir_pdbs_*.txt(N))
# [ "$pdb_logs" ] && rm $pdb_logs
java -jar $SRC/ProteinVolume_1.3.jar -het --radiusFileName $SRC/bondi.rad $1 |
    sed '1,6d' $1/OutputDir_pdbs_*.txt | sed $'s/   */\t/g' |
    mlr --tsv rename -g -r ' ,' then cut -x -f 'TimeTaken(ms)' then rename 'TotalVolume(A3),SolventExcludedVolume,VDWVolume,VanDerWaalsVolume' # stdout

