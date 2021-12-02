# GT
Prediction of reactivity for Glycosyl Transferase Superfamily 1 (GT1).

## REQUIREMENTS
- unix shell with zsh
- R
- pip
- miniconda for python 3

## OPTIONAL
- open source pymol (pymol.org) to run chemical feature generation pipeline. E.g. install with `conda install -c schrodinger pymol-bundle`

## CONFIG
- Configure some commands by running `./config.sh`.
- Add `src/` and `tools/degnlib/` to your `$PATH`.

## INSTALL
- `./install.sh` will install python packages to env "gt". Activate with `conda activate gt`.
