#!/usr/bin/env zsh
# make a command "git root" that gives the root folder of the repo.
git config alias.root 'rev-parse --show-toplevel'

# make sure you have miniconda installed
conda --version &> /dev/null || ./install_conda.sh

# Packages can be installing from requirements.txt with conda, then pip to set versions that were tested. 
# Easier install without version requirements using packages that are all in conda:
conda create --name GT --file requirements_conda.txt -c default -c conda-forge -c bioconda -c rdkit --yes

# Optionally install interactive python and support for neovim editing
conda install --name GT ipython
conda install --name GT -c conda-forge pynvim

# R packages
./install.R

# julia packages. Requires conda env with E3FP already installed, which is done above here.
./install.jl

