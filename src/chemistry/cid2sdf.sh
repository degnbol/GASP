#!/usr/bin/env zsh
# writes to multiple .sdf files named by the input CID
# cid2sdf.sh < cids.txt where cids.txt has a CID on each line

while read cid; do
	curl 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/'$cid'/record/SDF/?record_type=3d' > ${cid}.sdf 
done

