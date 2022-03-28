# GT
Prediction of reactivity for Glycosyl Transferase Superfamily 1 (GT1).

## REQUIREMENTS
- unix shell with zsh
- R
- pip
- miniconda for python 3, e.g. use `./install_conda.sh`

## OPTIONAL
- open source pymol (pymol.org) to run chemical feature generation pipeline. E.g. install with `conda install -c schrodinger pymol-bundle`

## CONFIG
- Configure `git root` commands by running `./config.sh`.
- Modify `$PATH` and `$PYTHONPATH` to easily access code, e.g. on zsh by running `./PATHS.sh >> ~/.zshrc`. Otherwise run `./PATHS.sh` and copy paste.

## INSTALL
- `./install.sh` will install python packages to env "gt". Activate with `conda activate gt`. It then calls `./install.R` which installs R packages.
- Optionally run `./results/generateResults.sh` to prepare some result files, e.g. feature data for training.

## RUN
- `encode_features.py --help` and `randomforest.py --help` for instructions.
- If using the same data from this work there is also `randomforest.sh` for convenience.

