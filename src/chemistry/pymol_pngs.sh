#!/usr/bin/env zsh
# Write .pngs with same filename as infiles in designated outdir.
# REQUIREMENT: pymol installed and accessible with command `pymol`
# USAGE: pymol_pngs.sh PDBs/ pngs/ [-H] [, width[, height[, dpi[, ray[, quiet]]]]]
# https://pymolwiki.org/index.php/Png
INDIR=$1
OUTDIR=$2
if [[ "$3" == "-H" ]]; then
    cmd0="hide (hydro); "
    EXTRA=${@:4}
else
    cmd0=""
    EXTRA=${@:3}
fi


# make sure dirs are supplied.
if [[ "${OUTDIR:e}" == "pdb" ]]; then
	echo "OUTDIR not supplied" 1>&2
	exit 1
fi

mkdir -p $OUTDIR

for fname in $INDIR/*.pdb; do
	name=${fname:t:r}
	echo $name
	
	cmd="${cmd0}png $OUTDIR/$name.png $EXTRA"
	pymol $fname -Q -cmd $cmd
done
