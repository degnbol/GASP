# GT
Prediction of reactivity for Glycosyl Transferase Superfamily 1 (GT1).

## REQUIREMENTS
- unix shell with zsh
- R
- pip
- miniconda for python 3, e.g. use `./install_conda.sh`

## INSTALL
- Modify `$PATH` and `$PYTHONPATH` to easily access code. Assuming zsh this is done by running `./PATHS.sh >> ~/.zshrc` or by copy-pasting the output of `./PATHS.sh` to somewhere in `~/.zshrc`.
- Run `./install.sh`. It will install python packages to environment "gt". Activate env with `conda activate gt`. `install.sh` also defines command `git root` and calls `./install.R` which installs R packages.
- Pymol in order to run chemical feature generation pipeline. E.g. install free open source pymol for [linux](https://pymolwiki.org/index.php/Linux_Install), [mac](https://pymolwiki.org/index.php/MAC_Install), or [windows](https://pymolwiki.org/index.php/Windows_Install). With homebrew on mac: `brew install brewsci/bio/pymol`
- Optionally run `./results/generateResults.sh` to prepare some result files, e.g. feature data for training.

## RUN
- `encode_features.py --help` and `randomforest.py --help` for instructions.
- If using the same data from this work there is also `randomforest.sh` for convenience.

