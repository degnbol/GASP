New chemical features can be generated given pubchem CIDs here.
See `acceptors_features.tsv.sh` for example of the original training set 
chemical features generation.
See `20220215_features.tsv.sh` for example of adding new chemicals.
The two things differ only by the fact that the 3D projection of chemicals will 
be different for a new set of chemicals so new ones needs to use the projection 
of the original set, otherwise all have to be run together.
