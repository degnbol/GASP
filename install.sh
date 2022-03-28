#!/usr/bin/env zsh
# make sure you have miniconda installed
conda --version &> /dev/null || ./install_conda.sh

# installing from requirements sets versions that was tested. There is some overlap in packages listed below.
conda create --name gt --file requirements.txt

# which packages are needed depends on what needs to be run.
# we will need muscle for multiple sequence alignment and hmmer to do HMM alignment on multiple sequence alignments
conda install -c bioconda muscle hmmer
# RDKit is for physiochemcial properties of acceptor molecules. https://www.rdkit.org/docs/Install.html
# It has to be a specific version so maybe have an isolated env for it
conda install -c rdkit rdkit
# 3D fingerprint
conda install -c conda-forge e3fp
# optionals for e3fp
conda install -c conda-forge h5py standardiser

# if we want to give attention to the balancing of positive and negative data points.
conda install imbalanced-learn
conda install cirpy # https://cirpy.readthedocs.io/en/latest/ conversion of CAS and other identifiers
conda install multiprocess # fork of multiprocessing that uses dill instead of pickle which makes it possible to async call local functions
