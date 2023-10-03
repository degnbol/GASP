# GASP
Glycosyltransferase Acceptor Specificity Prediction (GASP).
Pan-specific prediction of reactivity between any Glycosyltransferase 
Superfamily 1 (GT1) and chemical acceptor.

The install instructions here are thorough in order to allow for reproducing 
any of the work.
The main randomforest train and test script is `src/randomforest.py`.

## REQUIREMENTS
- Unix shell with zsh.
- R.
- Java runtime for some chemical feature generation.
  E.g. install with homebrew on Mac: `brew install adoptopenjdk8`
- `pip`.
- miniconda for python 3, e.g. use `./install_conda.sh`
- julia for
  - adding new chemicals to an already established E3FP MDS (chemical 
    features).
  - Feature selection (v1.8.5 used)

## INSTALL
- Modify `$PATH` and `$PYTHONPATH` to easily access code. Assuming zsh this is 
  done by running `./PATHS.sh >> ~/.zshrc` or by copy-pasting the output of 
  `./PATHS.sh` to somewhere in `~/.zshrc`.
- Run `./install.sh`. It will install python packages to environment "GT". 
  Activate env with `conda activate GT`.
  `install.sh` also defines command `git root` and calls `./install.R` which 
  installs R packages.
- [Miller](https://github.com/johnkerl/miller), e.g. with homebrew on mac: 
  `brew install miller`
- Pymol in order to run chemical feature generation pipeline. E.g. install free 
  open source pymol for [linux](https://pymolwiki.org/index.php/Linux_Install), 
  [mac](https://pymolwiki.org/index.php/MAC_Install), or 
  [windows](https://pymolwiki.org/index.php/Windows_Install). With homebrew on 
  mac: `brew install brewsci/bio/pymol`

## RUN
- Adding new acceptors: see `results/5-chemicalFeatures/`
- Adding new enzymes: see `results/6-unaligned/` and/or `results/7-align/`
- `encode_features.py --help` and `randomforest.py --help` for instructions.

